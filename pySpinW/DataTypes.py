import numpy as np
import dask.array as da
from .MatlabProxyObject import MatlabProxyObject
from .MatlabFunction import MatlabFunction

class DataTypes:

    def __init__(self, interface, pyMatlab, daskArray=False, chunks=None):
        """
        Data Converter to/from python and MATLAB
        :param matlab:
        """
        self._matlab = pyMatlab
        self._interface = interface
        self._dask = daskArray
        self._chunks = chunks
        # self.outNumpy = True
        # self.transpose = True

    def encode(self, data):

        # What is data?
        # 1) If it's a numpy array or a list, we convert. to matlab.double
        # 2) If it is a dict, then we enumerate values and encode them.
        # 3) If it's a tuple, it's a cell, which we enumerate. BUT, then we convert it into a list.
        # 4) If it is a double it's a double, if a integer, we encode to a double as well. MATLAB is tricky :-/

        if isinstance(data, (list, np.ndarray)):
            # Case 1)
            if isinstance(data, np.ndarray):
                if np.iscomplexobj(data):
                    data = self._interface.call('complex', (self.encode(data.real.tolist()), self.encode(data.imag.tolist())), nargout=1)
                else:
                    data = data.tolist()
                    data = self._matlab.double(data)
            else:
                data = self._matlab.double(data)
        elif isinstance(data, np.integer):
            # Case 4)
            data = float(data)
        elif isinstance(data, np.double):
            # Case 4)
            pass
        elif isinstance(data, MatlabFunction):
            data = data._fun
        elif isinstance(data, MatlabProxyObject):
            data = data.handle
        else:
            # Case 2, 3
            if isinstance(data, dict):
                # Case 2)
                for key, item in data.items():
                    data[key] = self.encode(item)
            elif isinstance(data, tuple):
                # Case 3)
                newdata = []
                for item in data:
                    newdata.append(self.encode(item))
                data = newdata

        # Unknown data i.e. text should pass through
        # TODO Make sure this works for more data cases...
        return data

    def decode(self, data):
        # This is where multiple results are passed back
        if isinstance(data, tuple):
            outData = []
            for thisData in data:
                outData.append(self.decode(thisData))
            return tuple(outData)
        # Decode the numeric data types. NOTE that we let the functions/methods slip through.
        if isinstance(data, list):
            # This is a cell return
            for key, item in enumerate(data):
                data[key] = self.decode(data[key])
            data = tuple(data)
        elif isinstance(data, self._matlab.double):
            data = self._arrayCreator(data)
        elif isinstance(data, self._matlab.int8):
            # TODO for all available data types
            data = self._arrayCreator(data, dtype=np.int)
        elif isinstance(data, str):
            if len(data) == 34:
                if data[0:2] == '!$':
                    try:
                        data = self.decode(self._interface.get_global(data[2:]))
                    except Exception as e:
                        print(e)
        elif isinstance(data, dict):
            for key, item in data.items():
                data[key] = self.decode(item)
        elif self._interface.callObj('isobject', data, nargout=1):
            data = MatlabProxyObject(self._interface, data, self)
        elif self._interface.call('isa', (data, 'function_handle'), nargout=1):
            data = MatlabFunction(self._interface, data, converter=self, parent=[])

        return data

    def _arrayCreator(self, data, dtype=None):
        if self._dask:
            return da.from_array(np.array(data), self._chunks)
        else:
            if dtype is not None:
                return np.array(data)
            else:
                return np.array(data, dtype=dtype)

    def _arrayReturner(self, data):
        if isinstance(data, np.ndarray):
            return data.tolist()
        elif isinstance(data, da.array):
            print('---- This is not advised -----')
            data =  np.array(data).tolist()