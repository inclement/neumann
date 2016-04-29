

with open('main.py') as fileh:
    main_string = fileh.read()

with open('neumann.kv') as fileh:
    kv_string = fileh.read()

kv_replacement = """
from kivy.lang import Builder
Builder.load_string('''
{}
''')
""".format(kv_string)

combined_string = main_string.replace('# kv goes here', kv_replacement.replace('\\n', '\\\n'))

with open('neumann_gui.py', 'w') as fileh:
    fileh.write(combined_string)
                                      
