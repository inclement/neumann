from kivy.uix.floatlayout import FloatLayout
from kivy.graphics import ClearColor, ClearBuffers, Rectangle, Canvas, Fbo
from kivy.properties import StringProperty, ListProperty, ObjectProperty
from kivy.clock import Clock
from kivy.lang import Builder

Builder.load_string('''
<ShaderWidget>:
    canvas:
        Color:
            rgba: 1, 1, 1, 1
        Rectangle:
            texture: self.fbo.texture if self.fbo else None
            tex_coords: [0, 0, 1, 0, 1, 1, 0, 1]
            pos: self.pos
            size: self.size
''')

class ShaderWidget(FloatLayout):

    # property to set the source code for fragment shader
    fs = StringProperty(None)

    # texture of the framebuffer
    texture = ObjectProperty(None)

    fbo_size = ListProperty([10, 10])

    fbo = ObjectProperty(None)

    def on_fbo_size(self, *args):
        self.fbo.size = self.fbo_size

    def __init__(self, **kwargs):
        self.canvas = Canvas()

        with self.canvas: 
            self.fbo = Fbo(size=self.fbo_size, use_parent_projection=True)

        with self.fbo:
            ClearColor(0, 0, 0, 0)
            ClearBuffers()
            Rectangle(size=(10000, 10000))

        super(ShaderWidget, self).__init__(**kwargs)

        self.update_glsl()
        #Clock.schedule_interval(self.update_glsl, 0)

        self.texture = self.fbo.texture

    def update_glsl(self, *largs):
        self.fbo['resolution'] = map(float, self.fbo_size)

    def on_fs(self, instance, value):
        #set the fragment shader to our source code
        shader = self.fbo.shader
        old_value = shader.fs
        shader.fs = value
        if not shader.success:
            shader.fs = old_value
            raise Exception('failed')
        self.update_glsl()

    #
    # now, if we have new widget to add,
    # add their graphics canvas to our Framebuffer, not the usual canvas.
    #

    def add_widget(self, widget):
        c = self.canvas
        self.canvas = self.fbo
        super(ShaderWidget, self).add_widget(widget)
        self.canvas = c

    def remove_widget(self, widget):
        c = self.canvas
        self.canvas = self.fbo
        super(ShaderWidget, self).remove_widget(widget)
        self.canvas = c

