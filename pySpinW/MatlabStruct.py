# from collections.abc import MutableMapping
#
# class MatlabStruct(MutableMapping):
#
#     def __init__(self, *args, interface=None, parent=None, converter=None, **kwargs):
#         self.parent = parent
#         self.interface = interface
#         self.converter = converter
#         self.__dict__.update(*args, **kwargs)
#
#     def __setitem__(self, key, value):
#         if self.interface is not None:
#             if self.interface.fieldnames(self.parent)
#         access = self.interface.substruct('.', key)
#         self.__dict__['handle'] = self.interface.subsasgn(self.handle, access, self.converter.encode(value))
#         self.__dict__[key] = value
#
#     def __getitem__(self, key):
#         return self.__dict__[key]
#
#     def __delitem__(self, key):
#         del self.__dict__[key]
#
#     def __iter__(self):
#         return iter(self.__dict__)
#
#     def __len__(self):
#         return len(self.__dict__)
#
#     def __str__(self):
#         '''returns simple dict representation of the mapping'''
#         return str(self.__dict__)
#     def __repr__(self):
#         '''echoes class, id, & reproducible representation in the REPL'''
#         return '{}, D({})'.format(super(D, self).__repr__(),
#                                   self.__dict__)