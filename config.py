#!/usr/bin/env python
'''

Part of parallyze
by Elizabeth Baird

this file originates in peasoup
by Camille Scott
@ www.github.com/camillescott/peasoup

Defines a simple config file parsing class.

'''

def coroutine(func):
    def start(*args,**kwargs):
        cr = func(*args,**kwargs)
        cr.next()
        return cr
    return start

class SimpleConfig(object):

    def __init__(self, key_types, fn):
        self.key_types = key_types
        self.fn = fn
        self.parse()
    
    def parse(self):
        
        parser = SimpleConfig.parse_protocol(
                key_types=self.key_types, target=self.parser_receiver())
        with open(self.fn, 'rb') as fp:
            for line in fp:
                for symbol in line.strip().split():
                    parser.send(symbol.strip())
            parser.send(None)

    @staticmethod
    @coroutine
    def parse_protocol(key_types=None, target=None):

        key = None
        value = []

        while True:
            symbol = (yield)
            if symbol in key_types or symbol == None:
                if value:
                    target.send((key,value))
                    value = []
                key = symbol
            elif key:
                value.append(symbol)
    
    @coroutine
    def parser_receiver(self):

        while True:
            key, value = (yield)
            value = [self.key_types[key](v) for v in value]
            if len(value) == 1:
                value = value[0]

            print '{k} : {v}'.format(k=key, v='\n'.join(value) \
                                    if type(value)==list else value)
            
            setattr(self, key, value)
            
            
