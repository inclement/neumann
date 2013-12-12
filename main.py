'''A speculative implementation of Neumann domain detection on the
gpu, via a kivy app.

'''

from kivy.app import App
from kivy.properties import (StringProperty, ListProperty, NumericProperty,
                             ObjectProperty, BooleanProperty, NumericProperty)
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.popup import Popup
from kivy.uix.scrollview import ScrollView
from kivy.clock import Clock
from kivy.animation import Animation

from shaderwidget import ShaderWidget

import random
from itertools import count
from functools import partial
from math import sqrt, pi

__version__ = '0.2'

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

shader_uniforms = '''
uniform vec2 resolution;
uniform float period;
uniform float max_intensity;
uniform float gradient_opacity;
uniform float phase_increment;
uniform float nearbys;
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

shader_bottom = '''

    return value;
}

vec2 grad(float pos_x, float pos_y, float dr) {
    float dfdx = (superposition_function(pos_x, pos_y) - superposition_function(pos_x + dr, pos_y)) / dr;
    float dfdy = (superposition_function(pos_x, pos_y) - superposition_function(pos_x, pos_y + dr)) / dr;
    return vec2(dfdx, dfdy);
}

mat2 rotmat(float angle) {
    return mat2(cos(angle), -1.0*sin(angle), sin(angle), cos(angle));
}

float mag(vec2 vector) {
    return sqrt(vector.x*vector.x + vector.y*vector.y);
}

float gradient(float pos_x, float pos_y, float dr)
{
    vec2 current_gradient = grad(pos_x, pos_y, dr);

    return atan(current_gradient.y, current_gradient.x) + 3.14159265;
}

/*
vec3 get_correlation_at(float pos_x, float pos_y, float dr) {
    vec2 local_grad = grad(pos_x, pos_y, dr);
    vec2 orth_dir = rotmat(3.14159/2.0) * local_grad;
    orth_dir = orth_dir / mag(orth_dir);

}
*/

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float dr = period / resolution.x;

    float value = superposition_function(pos_x, pos_y);

    vec3 intensity_col;
    vec3 gradient_col;
    vec3 correlation_col;
    intensity_col = intensity_to_colour(value / max_intensity);
    gradient_col = hsv2rgb( vec3(gradient(pos_x, pos_y, dr) / (2.0*3.1416), 1.0, 1.0));
    
    /*
    float corr = get_correlation_at(pos_x, pos_y, dr);
    float correlation_col = vec3(1.0-corr, 1.0-corr, 1.0-corr);
*/

    gl_FragColor = vec4(gradient_col.x*gradient_opacity + intensity_col.x*(1.0-gradient_opacity),
                        gradient_col.y*gradient_opacity + intensity_col.y*(1.0-gradient_opacity),
                        gradient_col.z*gradient_opacity + intensity_col.z*(1.0-gradient_opacity),
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
    plasma_shader = header + shader_uniforms + fileh.read()

with open('gradient.glsl') as fileh:
    gradient_shader = header + simple_shader_uniforms + fileh.read()

with open('tunnel_fly.glsl') as fileh:
    tunnel_fly_shader = header + shader_uniforms + fileh.read()

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

class AdvancedShader(ShaderWidget):

    shader_mid = StringProperty('')
    shader_uniforms = StringProperty('')
    fbo_texture = ObjectProperty(None, allownone=True)


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

    def replace_shader(self, *args):
        new_fs = (header + shader_uniforms + self.shader_parameters +
                  shader_top + self.shader_mid + shader_bottom)
        print 'new_fs is', new_fs
        self.fs = new_fs

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
        self.set_periodic_shader(scale=17, number=50)
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
    def on_gradient_opacity(self, *args):
        self.update_glsl()

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
        self.fbo['gradient_opacity'] = self.gradient_opacity
        self.fbo['phase_increment'] = float(self.phase_increment)
        self.fbo['nearbys'] = float(1.0)

    def view_wavevectors(self, *args):
        WvPopup(content=WvPopupContent(wavevectors=self.wavevectors)).open()

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
