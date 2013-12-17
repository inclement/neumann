'''A speculative implementation of Neumann domain detection on the
gpu, via a kivy app.

'''

from kivy.app import App
from kivy.properties import (StringProperty, ListProperty, NumericProperty,
                             ObjectProperty, BooleanProperty, NumericProperty)
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.popup import Popup
from kivy.uix.scrollview import ScrollView
from kivy.uix.widget import Widget
from kivy.clock import Clock
from kivy.animation import Animation
from kivy.graphics import BindTexture

from shaderwidget import ShaderWidget

import random
from itertools import count
from functools import partial
from math import sqrt, pi

__version__ = '0.3'

header = '''
#ifdef GL_ES
precision highp float;
#endif

/* Outputs from the vertex shader */
varying vec4 frag_color;
varying vec2 tex_coord0;

/* uniform texture samplers */
uniform sampler2D texture0;
'''

simple_shader_uniforms = '''
uniform vec2 resolution;
'''

universal_shader_uniforms = '''
uniform vec2 resolution;
'''
secondary_shader_uniforms = '''
uniform sampler2D input_texture;
uniform float fbo_jump;
uniform int while_cutoff;
'''

critical_finder_shader = (header + universal_shader_uniforms +
                          secondary_shader_uniforms) + \
'''
vec2 grad(float pos_x, float pos_y, float dr) {
    float dfdx = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x + dr, pos_y)) / dr;
    float dfdy = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x, pos_y + dr)) / dr;
    return vec2(dfdx, dfdy);
}

float gradient(float pos_x, float pos_y, float dr)
{
    vec2 current_gradient = grad(pos_x, pos_y, dr);

    return atan(current_gradient.y, current_gradient.x) + 3.14159265;
}
void main (void) {
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;
    float pos_x = x / resolution.x;
    float pos_y = y / resolution.y;
}
'''

correlation_shader_uniforms = '''
uniform float correlation_width;
uniform float orth_jump;
uniform float correlation_cutoff;
'''

neumann_shader_uniforms = '''
uniform float period;
uniform float max_intensity;
uniform float phase_increment;
'''

critical_shader_uniforms = '''
uniform float critical_jump;
'''

inversion_shader = header + universal_shader_uniforms + '''
uniform sampler2D input_texture;

void main(void) {
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;
    float pos_x = x / resolution.x;
    float pos_y = y / resolution.y;
    vec4 cin = texture2D(input_texture, vec2(pos_x, pos_y));
    gl_FragColor = vec4(1.0-cin.x, 1.0-cin.y, 1.0-cin.z, 1.0);
}
'''

shader_top = '''
vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 intensity_to_colour(float intensity)
{
    if (intensity > 0.0) {
        return vec3(1.0, 1.0 - intensity, 1.0 - intensity);
        }
    else {
        return vec3(1.0 - abs(intensity), 1.0 - abs(intensity), 1.0);
        }
}


float superposition_function(float pos_x, float pos_y)
{
    float value = 0.0;

'''

correlation_shader_bottom = '''
    return value;
}


mat2 rotmat(float angle) {
    return mat2(cos(angle), -1.0*sin(angle), sin(angle), cos(angle));
}

float mag(vec2 vector) {
    return sqrt(vector.x*vector.x + vector.y*vector.y);
}

vec2 grad(float pos_x, float pos_y, float dr) {
    float dfdx = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x + dr, pos_y)) / dr;
    float dfdy = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x, pos_y + dr)) / dr;
    return vec2(dfdx, dfdy);
}

float gradient(float pos_x, float pos_y, float dr)
{
    vec2 current_gradient = grad(pos_x, pos_y, dr);

    return atan(current_gradient.y, current_gradient.x) + 3.14159265;
}

float get_correlation_at(float pos_x, float pos_y, float dr) {
    vec2 local_grad = grad(pos_x, pos_y, dr);
    vec2 orth_dir = rotmat(3.14159/2.0) * local_grad;
    orth_dir = orth_dir / mag(orth_dir);

    vec2 orth_grad_1 = grad(pos_x - orth_dir.x * orth_jump,
                            pos_y - orth_dir.y * orth_jump, dr);
    vec2 orth_grad_2 = grad(pos_x + orth_dir.x * orth_jump,
                            pos_y + orth_dir.y * orth_jump, dr);

    float orth_angle_1 = atan(orth_grad_1.y, orth_grad_1.x);
    float orth_angle_2 = atan(orth_grad_2.y, orth_grad_2.y);
    float angle_grad = atan(local_grad.y, local_grad.x);

    float relative_angle_1 = angle_grad - orth_angle_1;
    if (relative_angle_1 > 3.14159) {
        relative_angle_1 -= 3.14159;
    }
    float relative_angle_2 = angle_grad - orth_angle_2;
    if (relative_angle_2 > 3.14159) {
        relative_angle_2 -= 3.14159;
    }

    float angle_sum = relative_angle_1 + relative_angle_2;
    return exp(-1.0 * (angle_sum * angle_sum) /
               (correlation_width*correlation_width));
}

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float dr = period / resolution.x;

    float corr = get_correlation_at(pos_x, pos_y, dr);
    float corrcol = 1.0 - corr;
    float cutoff_colour;

    if (corrcol < correlation_cutoff) {
        cutoff_colour = 0.0;
    } else {
        cutoff_colour = 1.0;
    }

    gl_FragColor = vec4(cutoff_colour, cutoff_colour,
                        cutoff_colour, 1.0 - cutoff_colour);
}
'''

gradient_shader_bottom = '''
    return value;
}

vec2 grad(float pos_x, float pos_y, float dr) {
    float dfdx = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x + dr, pos_y)) / dr;
    float dfdy = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x, pos_y + dr)) / dr;
    return vec2(dfdx, dfdy);
}

float gradient(float pos_x, float pos_y, float dr)
{
    vec2 current_gradient = grad(pos_x, pos_y, dr);

    return atan(current_gradient.y, current_gradient.x) + 3.14159265;
}

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float dr = period / resolution.x;

    vec3 gradient_col;
    gradient_col = hsv2rgb( vec3(gradient(pos_x, pos_y, dr) /
                            (2.0*3.1416), 1.0, 1.0));

    gl_FragColor = vec4(gradient_col.x,
                        gradient_col.y,
                        gradient_col.z,
                        1.0);
}
'''

critical_shader_bottom = '''
    return value;
}

vec2 grad(float pos_x, float pos_y, float dr) {
    float dfdx = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x + dr, pos_y)) / dr;
    float dfdy = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x, pos_y + dr)) / dr;
    return vec2(dfdx, dfdy);
}

float gradient(float pos_x, float pos_y, float dr)
{
    vec2 current_gradient = grad(pos_x, pos_y, dr);

    return atan(current_gradient.y, current_gradient.x) + 3.14159265;
}

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float dr = period / resolution.x;

    float angle_right = gradient(pos_x + critical_jump, pos_y, dr);
    float angle_down = gradient(pos_x, pos_y - critical_jump, dr);
    float angle_left = gradient(pos_x - critical_jump, pos_y, dr);
    float angle_up = gradient(pos_x, pos_y + critical_jump, dr);

    float d1 = angle_down - angle_right;
    float d2 = angle_left - angle_down;
    float d3 = angle_up - angle_left;
    float d4 = angle_right - angle_up;

    float pi = 3.14159;
    float twopi = 2.0*pi;

    if (d1 > pi) {
        d1 -= twopi;
    }
    if (d1 < -1.0 * pi) {
        d1 += twopi;
    }
    if (d2 > pi) {
        d2 -= twopi;
    }
    if (d2 < -1.0 * pi) {
        d2 += twopi;
    }
    if (d3 > pi) {
        d3 -= twopi;
    }
    if (d3 < -1.0 * pi) {
        d3 += twopi;
    }
    if (d4 > pi) {
        d4 -= twopi;
    }
    if (d4 < -1.0 * pi) {
        d4 += twopi;
    }

    float value = superposition_function(pos_x, pos_y);

    vec4 output_colour;
    if (all(bvec4(d1 > 0.0, d2 > 0.0, d3 > 0.0, d4 > 0.0))) {
        output_colour = vec4(1.0, 0.7, 0.0, 1.0);
    } else {
        if (all(bvec4(d1 < 0.0, d2 < 0.0, d3 < 0.0, d4 < 0.0))) {
            if (value < 0.0) {
                output_colour = vec4(0.0, 0.0, 0.8, 1.0);
            } else {
                output_colour = vec4(0.8, 0.0, 0.0, 1.0);
            }
        } else {
            output_colour = vec4(0.0, 0.0, 0.0, 0.0);

        }
    }

    gl_FragColor = output_colour;
}
'''

line_shader_bottom = '''
    return value;
}

vec2 grad(float pos_x, float pos_y, float dr) {
    float dfdx = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x + dr, pos_y)) / dr;
    float dfdy = (superposition_function(pos_x, pos_y) -
                  superposition_function(pos_x, pos_y + dr)) / dr;
    return vec2(dfdx, dfdy);
}

float gradient(float pos_x, float pos_y, float dr)
{
    vec2 current_gradient = grad(pos_x, pos_y, dr);

    return atan(current_gradient.y, current_gradient.x) + 3.14159265;
}

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float frac_x = x / resolution.x;
    float frac_y = y / resolution.y;

    float value = superposition_function(pos_x, pos_y);
    float dr = period / resolution.x;

    vec4 cin = texture2D(input_texture, vec2(frac_x, frac_y));

    vec4 current_colour;
    float cur_frac_x = frac_x;
    float cur_frac_y = frac_y;
    float cur_gradient;
    current_colour = texture2D(input_texture, vec2(frac_x, frac_y));
    int num_steps = 0;

    while (all(bvec2(current_colour.w < 0.9, num_steps < while_cutoff))) {
        cur_gradient = gradient(cur_frac_x * period, cur_frac_y * period, dr);
        cur_frac_x += 1.0*fbo_jump * cos(cur_gradient);
        cur_frac_y += 1.0*fbo_jump * sin(cur_gradient);
        current_colour = texture2D(input_texture, vec2(cur_frac_x, cur_frac_y));
        num_steps += 1;
    }

    vec4 output_colour;

    if (current_colour.x > 0.9) {
        output_colour = vec4(0.5, 0.0, 0.5, 1.0);
    } else {
        cur_frac_x = frac_x;
        cur_frac_y = frac_y;
        current_colour = texture2D(input_texture, vec2(frac_x, frac_y));
        num_steps = 0;
        while (all(bvec2(current_colour.w < 0.9, num_steps < while_cutoff))) {
            cur_gradient = gradient(cur_frac_x * period,
                                    cur_frac_y * period, dr);
            cur_frac_x -= 1.0*fbo_jump * cos(cur_gradient);
            cur_frac_y -= 1.0*fbo_jump * sin(cur_gradient);
            current_colour = texture2D(input_texture,
                                       vec2(cur_frac_x, cur_frac_y));
            num_steps += 1;
        }

        if (current_colour.x > 0.9) {
            output_colour = vec4(0.5, 0.0, 0.5, 1.0);
        } else {
            output_colour = vec4(0.0, 0.0, 0.0, 0.0);
        }
    }

    gl_FragColor = output_colour;
}
'''

intensity_shader_bottom = '''
    return value;
}

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float dr = period / resolution.x;

    float value = superposition_function(pos_x, pos_y);

    vec3 intensity_col;
    intensity_col = intensity_to_colour(value / max_intensity);

    gl_FragColor = vec4(intensity_col.x,
                        intensity_col.y,
                        intensity_col.z,
                        1.0);

}
'''

trivial_shader = header + '''
void main(void)
{
gl_FragColor = vec4(0.5, 0.2, 0.2, 1.0);
}
'''

with open('plasma.glsl') as fileh:
    plasma_shader = header + universal_shader_uniforms + fileh.read()

with open('gradient.glsl') as fileh:
    gradient_shader = header + simple_shader_uniforms + fileh.read()

with open('tunnel_fly.glsl') as fileh:
    tunnel_fly_shader = header + universal_shader_uniforms + fileh.read()

def duofactors(k):
    outs = []
    for i in range(int(sqrt(k))+2):
        for j in range(i, int(sqrt(k))+2):
            if (i**2 + j**2) == k:
                outs.append((i, j))
    return outs

def get_periodic_wavevectors(number=50, scale=5, seed=0):
    seed = random.randint(0, 1000000)
    generator = random.Random()
    generator.seed(seed)

    possible_wvs = duofactors(scale)
    possible_signs = [[1, 1], [1, -1], [-1, 1], [-1, -1]]

    amps = [generator.gauss(0, 1) for i in range(number)]
    phases = [2*pi*generator.random() for i in range(number)]
    wvs = [] #n.zeros((number, 2), dtype=n.float64)

    for i in range(number):
        wv = map(float, generator.choice(possible_wvs))
        generator.shuffle(wv)
        random_sign = generator.choice(possible_signs)
        wv[0] *= random_sign[0]
        wv[1] *= random_sign[1]
        wv = [2*pi/sqrt(scale/2.) * entry for entry in wv]
        wvs.append(wv)

    return list(wvs), phases, amps


class VSeparator(Widget):
    pass


class AdvancedShader(ShaderWidget):

    shader_mid = StringProperty('')
    shader_parameters = StringProperty('')
    fbo_texture = ObjectProperty(None, allownone=True)

    def __init__(self, *args, **kwargs):
        self.register_event_type('on_update_glsl')
        super(AdvancedShader, self).__init__(*args, **kwargs)

    def on_fs(self, *args):
        super(AdvancedShader, self).on_fs(*args)
        self.fbo_texture = self.fbo.texture

    def on_fbo_size(self, *args):
        super(AdvancedShader, self).on_fbo_size(*args)
        self.fbo.size = self.fbo_size
        self.fbo_texture = self.fbo.texture
        self.update_glsl()

    def update_glsl(self, *args):
        super(AdvancedShader, self).update_glsl(*args)

        self.fbo['resolution'] = map(float, self.fbo_size)
        self.dispatch('on_update_glsl')

    def on_update_glsl(self, *args):
        pass

    def replace_shader(self, *args):
        new_fs = (header + universal_shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + shader_bottom)
        self.fs = new_fs


class SecondaryShader(AdvancedShader):
    input_texture = ObjectProperty(None, allownone=True)
    parent_shader = ObjectProperty(None, allownone=True)
    search_distance = NumericProperty(0.02)
    def __init__(self, *args, **kwargs):
        super(SecondaryShader, self).__init__(*args, **kwargs)
        with self.fbo.before:
            self.input_binding = BindTexture(texture=self.input_texture,
                                             index=1)
        self.fbo['input_texture'] = 1
        self.update_binding()
        Clock.schedule_interval(self.force_update, 1/60.)
        Clock.schedule_interval(self.update_glsl, 1.)

    def on_fbo_size(self, *args):
        super(SecondaryShader, self).on_fbo_size(*args)
        self.fbo['fbo_size'] = map(float, self.fbo_size)

    def on_parent_shader(self, *args):
        parent_shader = self.parent_shader
        parent_shader.bind(on_update_glsl=self.force_update)

    def update_binding(self, *args):
        self.input_binding.texture = self.input_texture

    def on_input_texture(self, *args):
        self.update_binding()

    def force_update(self, *args):
        self.input_binding.force_update()

    def update_glsl(self, *args):
        super(SecondaryShader, self).update_glsl(*args)
        self.fbo['fbo_size'] = map(float, self.fbo_size)
        self.fbo['fbo_jump'] = float(1./self.fbo_size[0])


class NeumannShader(AdvancedShader):
    wavevectors = ListProperty([])
    phases = ListProperty([])
    period = NumericProperty(1.0)
    scale = NumericProperty(0)
    number = NumericProperty(0)
    downscale = NumericProperty(0)
    gradient_opacity = NumericProperty(0.0)
    shader_parameters = StringProperty('')
    phase_increment = NumericProperty(0.0)
    animation = ObjectProperty(None, allownone=True)

    def __init__(self, *args, **kwargs):
        super(NeumannShader, self).__init__(*args, **kwargs)
        self.set_periodic_shader(scale=17, number=25)

    def on_phase_increment(self, *args):
        self.fbo['phase_increment'] = float(self.phase_increment)

    def animate_phase_increment(self, time):
        if self.animation is not None:
            self.stop_animation()
            return
        anim = Animation(phase_increment=2*pi,
                         t='linear',
                         duration=(1 - self.phase_increment/(2*pi)) * time)
        anim.bind(on_complete=partial(self.repeat_animation, time))
        self.animation = anim
        anim.start(self)

    def repeat_animation(self, time, *args):
        self.phase_increment = 0.
        self.animation = None
        self.animate_phase_increment(time)

    def stop_animation(self, *args):
        if self.animation is not None:
            self.animation.cancel(self)
            self.animation = None

    def replace_shader(self, *args):
        new_fs = (header + universal_shader_uniforms + neumann_shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + intensity_shader_bottom)
        self.fs = new_fs

    def set_periodic_shader(self, scale=5, number=10, downscale=2):
        try:
            scale = int(scale)
        except ValueError:
            return  # If user entered a wrong value
        if not duofactors(scale):
            return False

        self.scale = scale
        self.number = number
        self.downscale = downscale

        length = int(100/float(downscale) * float(scale)/5.)
        #self.fbo_size = (length, length)
        # print 'length is', length
        # if self.parent:
            # slider = self.parent.resolution_slider
            # print 'slider is', slider
            # if slider.max < length:
            #     slider.max = length+1
            # self.parent.resolution_slider.value = length

        self.period = float(sqrt(scale/2.))

        new_wavevectors, new_phases, new_amplitudes = get_periodic_wavevectors(number, scale) 
        self.wavevectors = zip(map(list, new_wavevectors), map(float, new_phases), new_amplitudes)

    def on_wavevectors(self, *args):
        shader_mid = ''
        shader_parameters = ''
        i = 0
        for wv in self.wavevectors:
            current_uniform = 'k{}'.format(i)
            current_phase = 'phase{}'.format(i)
            current_amp = 'amp{}'.format(i)
            shader_parameters += ('''
            uniform vec2 {};
            uniform float {};
            uniform float {};
            ''').format(current_uniform, current_phase, current_amp)
            shader_mid += ('''
            value += {amp}*sin({cu}.x * pos_x + {cu}.y * pos_y + {ph} + phase_increment);
            ''').format(amp=current_amp, cu=current_uniform, ph=current_phase)
            i += 1
        self.shader_parameters = shader_parameters
        self.shader_mid = shader_mid
        self.replace_shader()
        self.update_glsl()

    def update_glsl(self, *args):
        super(NeumannShader, self).update_glsl(*args)
        for number, param in zip(count(), self.wavevectors):
            wv, phase, amp = param
            current_wv = 'k{}'.format(number)
            current_phase = 'phase{}'.format(number)
            current_amp = 'amp{}'.format(number)
            self.fbo[current_wv] = [float(wv[0]), float(wv[1])]
            self.fbo[current_phase] = phase
            self.fbo[current_amp] = amp
            number += 1
        self.fbo['period'] = float(self.period)
        self.fbo['max_intensity'] = float(2*sqrt(max(1.0, float(len(self.wavevectors)))))
        self.fbo['phase_increment'] = float(self.phase_increment)

    def view_wavevectors(self, *args):
        WvPopup(content=WvPopupContent(wavevectors=self.wavevectors)).open()


class LineDetectionShader(NeumannShader):
    input_texture = ObjectProperty(None, allownone=True)
    parent_shader = ObjectProperty(None, allownone=True)
    search_distance = NumericProperty(0.02)
    while_cutoff = NumericProperty(1)

    def __init__(self, *args, **kwargs):
        super(LineDetectionShader, self).__init__(*args, **kwargs)
        with self.fbo.before:
            self.input_binding = BindTexture(texture=self.input_texture,
                                             index=1)
        self.fbo['input_texture'] = 1
        self.update_binding()
        Clock.schedule_interval(self.force_update, 1/60.)

    def on_while_cutoff(self, *args):
        self.fbo['while_cutoff'] = int(self.while_cutoff)

    def on_fbo_size(self, *args):
        super(LineDetectionShader, self).on_fbo_size(*args)
        self.fbo['fbo_size'] = map(float, self.fbo_size)
        self.fbo['fbo_jump'] = float(1./self.fbo_size[0])

    def on_parent_shader(self, *args):
        parent_shader = self.parent_shader
        parent_shader.bind(on_update_glsl=self.force_update)

    def update_binding(self, *args):
        self.input_binding.texture = self.input_texture

    def on_input_texture(self, *args):
        self.update_binding()

    def force_update(self, *args):
        self.input_binding.force_update()

    def update_glsl(self, *args):
        super(LineDetectionShader, self).update_glsl(*args)
        self.fbo['fbo_size'] = map(float, self.fbo_size)
        self.fbo['fbo_jump'] = float(1./self.fbo_size[0])
        self.fbo['while_cutoff'] = int(self.while_cutoff)

    def replace_shader(self, *args):
        new_fs = (header + universal_shader_uniforms + secondary_shader_uniforms +
                  neumann_shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + line_shader_bottom)
        print 'new_fs is'
        print new_fs
        self.fs = new_fs


class GradientShader(NeumannShader):
    def __init__(self, *args, **kwargs):
        super(GradientShader, self).__init__(*args, **kwargs)

    def replace_shader(self, *args):
        new_fs = (header + universal_shader_uniforms + neumann_shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + gradient_shader_bottom)
        self.fs = new_fs


class CriticalShader(GradientShader):
    critical_jump = NumericProperty(0.1)

    def __init__(self, *args, **kwargs):
        super(GradientShader, self).__init__(*args, **kwargs)

    def on_critical_jump(self, *args):
        self.fbo['critical_jump'] = float(self.critical_jump)

    def replace_shader(self, *args):
        new_fs = (header + universal_shader_uniforms + neumann_shader_uniforms + critical_shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + critical_shader_bottom)
        self.fs = new_fs

    def update_glsl(self, *args):
        super(CriticalShader, self).update_glsl(*args)
        self.fbo['critical_jump'] = float(self.critical_jump)


class CorrelationShader(NeumannShader):
    orth_jump = NumericProperty(0.001)
    correlation_width = NumericProperty(0.005)
    correlation_cutoff = NumericProperty()
    def __init__(self, *args, **kwargs):
        super(CorrelationShader, self).__init__(*args, **kwargs)

    def on_correlation_width(self, *args):
        self.fbo['correlation_width'] = float(self.correlation_width)

    def on_orth_jump(self, *args):
        self.fbo['orth_jump'] = float(self.orth_jump)

    def on_correlation_cutoff(self, *args):
        self.fbo['correlation_cutoff'] = float(self.correlation_cutoff)

    def replace_shader(self, *args):
        new_fs = (header + universal_shader_uniforms + correlation_shader_uniforms +
                  neumann_shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + correlation_shader_bottom)
        self.fs = new_fs

    def update_glsl(self, *args):
        super(CorrelationShader, self).update_glsl(*args)
        self.fbo['orth_jump'] = float(self.orth_jump)
        self.fbo['correlation_width'] = float(self.correlation_width)
        self.fbo['correlation_cutoff'] = float(self.correlation_cutoff)


class WvPopup(Popup):
    wavevectors = ListProperty([])


class WvPopupContent(BoxLayout):
    wavevectors = ListProperty([])


class Interface(BoxLayout):
    shader_display = ObjectProperty()
    resolution_slider = ObjectProperty()
    def __init__(self, *args, **kwargs):
        super(Interface, self).__init__(*args, **kwargs)


class NeumannApp(App):
    def build(self):
        return Interface()

    def on_pause(self):
        return True

if __name__ == "__main__":
    NeumannApp().run()
