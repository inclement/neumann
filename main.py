'''A speculative implementation of Neumann domain detection on the
gpu, via a kivy app.

'''

from kivy.app import App
from kivy.properties import StringProperty, ListProperty, NumericProperty, ObjectProperty
from kivy.uix.boxlayout import BoxLayout
from kivy.clock import Clock

from shaderwidget import ShaderWidget

import random
import numpy as n
from itertools import count

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
shader_uniforms = '''
uniform vec2 resolution;
uniform float period;
uniform float max_intensity;
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

void main(void)
{
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;

    float pos_x = x / resolution.x * period;
    float pos_y = y / resolution.y * period;

    float value = 0.0;

'''

shader_bottom = '''
   vec3 rgbcol = intensity_to_colour(value / max_intensity);

   gl_FragColor = vec4(rgbcol.x, rgbcol.y, rgbcol.z, 1.0);
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
    gradient_shader = header + shader_uniforms + fileh.read()

with open('tunnel_fly.glsl') as fileh:
    tunnel_fly_shader = header + shader_uniforms + fileh.read()

def duofactors(k):
    outs = []
    for i in range(int(n.sqrt(k))+2):
        for j in range(i, int(n.sqrt(k))+2):
            if (i**2 + j**2) == k:
                outs.append((i, j))
    return outs

def get_periodic_wavevectors(number=50, scale=5, seed=0):
    seed = random.randint(0, 1000000)
    generator = n.random.RandomState()
    generator.seed(seed)

    possible_wvs = duofactors(scale)
    possible_signs = map(n.array, [[1, 1], [1, -1], [-1, 1], [-1, -1]])

    amps = generator.normal(size=number)
    phases = 2*n.pi*generator.rand(number)
    wvs = n.zeros((number, 2), dtype=n.float64)

    for i in range(number):
        wv = n.array(possible_wvs[generator.randint(len(possible_wvs))],
                     dtype=n.float64)
        generator.shuffle(wv)
        wv *= possible_signs[generator.randint(len(possible_signs))]
        wv *= 2*n.pi/n.sqrt(scale/2.)
        wvs[i] = wv

    return list(wvs)


class NeumannShader(ShaderWidget):
    wavevectors = ListProperty([])

    shader_mid = StringProperty('')
    shader_uniforms = StringProperty('')
    period = NumericProperty(1.0)
    fbo_texture = ObjectProperty(None, allownone=True)

    scale = NumericProperty(0)
    number = NumericProperty(0)
    downscale = NumericProperty(0)

    def __init__(self, *args, **kwargs):
        super(NeumannShader, self).__init__(*args, **kwargs)
        self.set_periodic_shader(scale=65, number=50)

    def on_fs(self, *args):
        super(NeumannShader, self).on_fs(*args)
        # print 'new fs is:'
        # print self.fs

    def on_fbo_size(self, *args):
        super(NeumannShader, self).on_fbo_size(*args)
        self.fbo.size = self.fbo_size
        self.fbo_texture = self.fbo.texture
        self.update_glsl()

    def set_periodic_shader(self, scale=5, number=10, downscale=2):
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

        self.period = float(n.sqrt(scale/2.))

        new_wavevectors = get_periodic_wavevectors(number, scale) 
        print 'new_wavevectors are', new_wavevectors
        self.wavevectors = map(list, new_wavevectors)
        print 'set them.'

    def on_wavevectors(self, *args):
        shader_mid = ''
        shader_uniforms = ''
        i = 0
        for wv in self.wavevectors:
            current_uniform = 'k{}'.format(i)
            shader_uniforms += ('''
            uniform vec2 {};
            ''').format(current_uniform)
            shader_mid += ('''
            value += sin({cu}.x * pos_x + {cu}.y * pos_y);
            ''').format(cu=current_uniform)
            i += 1
        self.shader_uniforms = shader_uniforms
        self.shader_mid = shader_mid
        self.replace_shader()
        self.update_glsl()

    def update_glsl(self, *args):
        super(NeumannShader, self).update_glsl(*args)
        print 'self.fs is', self.fs

        self.fbo['resolution'] = map(float, self.fbo_size)
        print 'set resolution to', map(float, self.fbo_size)

        for number, wv in zip(count(), self.wavevectors):
            current_uniform = 'k{}'.format(number)
            self.fbo[current_uniform] = [float(wv[0]), float(wv[1])]
            number += 1
        self.fbo['period'] = float(self.period)
        print 'set period to', float(self.period)

        self.fbo['max_intensity'] = float(n.sqrt(max(1.0, float(len(self.wavevectors)))))
        print 'set max_intensity to', max(1.0, float(len(self.wavevectors))) 

    def replace_shader(self, *args):
        self.fs = (header + shader_uniforms + self.shader_uniforms +
                   shader_top + self.shader_mid + shader_bottom)

class Interface(BoxLayout):
    shader_display = ObjectProperty()
    resolution_slider = ObjectProperty()
    def __init__(self, *args, **kwargs):
        super(Interface, self).__init__(*args, **kwargs)

class NeumannApp(App):
    def build(self):
        return Interface()

if __name__ == "__main__":
    NeumannApp().run()
