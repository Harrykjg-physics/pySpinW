import re
from io import StringIO

from .MatlabFunction import MatlabFunction

class MatlabProxyObject:
    """A Proxy for an object that exists in Matlab.

    All property accesses and function calls are executed on the
    Matlab object in Matlab.

    Auto populates methods and properties and can be called ala matlab/python

    """

    def __init__(self, interface, handle, converter):
        """
        Create a non numeric object of class handle (an object from a non-numeric class).
        :param interface: The callable MATLAB interface (where we run functions)
        :param handle: The matlabObject which represents a class object
        :param converter: The converter class between MATLAB/python
        """
        self.__dict__['handle'] = handle
        self.__dict__['interface'] = interface
        self.__dict__['converter'] = converter

        for attribute in self._getAttributeNames():
            self.__dict__[attribute] = self.__getattr__(attribute)
        for method in self._getMethodNames():
            super(MatlabProxyObject, self).__setattr__(method,
                                                       MatlabFunction(self.interface, method,
                                                                      converter=self.converter, parent=self.handle,
                                                                      caller=self))

    def _getAttributeNames(self):
        """
        Gets attributes from a MATLAB object
        :return: list of attribute names
        """
        return self.interface.callObj('fieldnames', self.handle)

    def _getMethodNames(self):
        """
        Gets methods from a MATLAB object
        :return: list of method names
        """
        return self.interface.callObj('methods', self.handle)

    def __getattr__(self, name):
        """Retrieve a value or function from the object.

        Properties are returned as native Python objects or
        :class:`MatlabProxyObject` objects.

        Functions are returned as :class:`MatlabFunction` objects.

        """
        interface = self.interface
        # if it's a property, just retrieve it
        if name in interface.callObj('properties', self.handle, nargout=1):
            return interface.call('subsref', (self.handle, interface.call('substruct', ('.', name))))
        # if it's a method, wrap it in a functor
        elif name in interface.callObj('methods', self.handle, nargout=1):
            class matlab_method:
                def __call__(_self, *args, nargout=-1, **kwargs):
                    # serialize keyword arguments:
                    args += sum(kwargs.items(), ())
                    return getattr(interface, name)(self, *args, nargout=nargout)

                # only fetch documentation when it is actually needed:
                @property
                def __doc__(_self):
                    classname = getattr(interface, 'class')(self)
                    # TODO this is probably broken :-/
                    return interface.help('{0}.{1}'.format(classname, name), nargout=1)

            return matlab_method()
        else:
            raise AttributeError

    def __setattr__(self, name, value):
        self.__class__[name] = value
        access = self.handle, self.interface.call('substruct', ('.', name))
        self.interface.call('subsasgn', (self, access, value))

    def __repr__(self):
        # getclass = self.interface.str2func('class')
        return "<proxy for Matlab {} object>".format(self.interface.callObj('class', self.handle))

    def __str__(self):
        # remove pseudo-html tags from Matlab output
        # html_str = self.interface.call2('str2func', "@(x) evalc('disp(x)')")
        html_str = self.interface.callObj(self.interface.callObj('str2func', "@(x) evalc('disp(x)')"), self.handle)
        return re.sub('</?a[^>]*>', '', html_str)

    @property
    def __doc__(self):
        out = StringIO()
        return self.interface.help(self.handle, nargout=1, stdout=out)

    def updateProxy(self):
        """
        Perform a update on an objects fields. Useful for when dealing with handle classes.
        :return: None
        """
        # We assume methods can't change
        for attribute in self._getAttributeNames():
            self.__dict__[attribute] = self.__getattr__(attribute)

    def updateObj(self):
        """
        When you change an attributes value, the corresponding value should be changed in self.handle
        :return: None
        """
        raise NotImplementedError
