import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

class flowmap( dict ): pass
def flowmap_rep(dumper, data):
    return dumper.represent_mapping( u'tag:yaml.org,2002:map', data, flow_style=True)

class blockseqtrue( list ): pass
def blockseqtrue_rep(dumper, data):
    return dumper.represent_sequence( u'tag:yaml.org,2002:seq', data, flow_style=True )

yaml.add_representer(blockseqtrue, blockseqtrue_rep)
yaml.add_representer(flowmap, flowmap_rep)

class MyDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()

def FormatSettings_main(data):
    
    if "condensation-rate" in data['planet']['water'].keys():
        data['planet']['water']['condensation-rate'] = flowmap(data['planet']['water']['condensation-rate'])
    
    if 'particles' in data:
        for i in range(len(data['particles'])):
            if "condensation-rate" in data['particles'][i]:
                data['particles'][i]["condensation-rate"] = \
                flowmap(data['particles'][i]["condensation-rate"])

    for i in range(len(data['boundary-conditions'])):
        if "lower-boundary" in data['boundary-conditions'][i]:
            order = ['type','vdep','mix','flux','height']
            copy = data['boundary-conditions'][i]['lower-boundary'].copy()
            data['boundary-conditions'][i]['lower-boundary'].clear()
            for key in order:
                if key in copy.keys():
                    data['boundary-conditions'][i]['lower-boundary'][key] = copy[key]

            data['boundary-conditions'][i]['lower-boundary'] = flowmap(data['boundary-conditions'][i]['lower-boundary'])

            order = ['type','veff','flux']
            copy = data['boundary-conditions'][i]['upper-boundary'].copy()
            data['boundary-conditions'][i]['upper-boundary'].clear()
            for key in order:
                if key in copy.keys():
                    data['boundary-conditions'][i]['upper-boundary'][key] = copy[key]

            data['boundary-conditions'][i]['upper-boundary'] = flowmap(data['boundary-conditions'][i]['upper-boundary'])
        
    return data
    
    