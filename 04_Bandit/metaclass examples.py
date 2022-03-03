# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:bandit_38]
#     language: python
#     name: conda-env-bandit_38-py
# ---

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %% [markdown]
# ## Metaclass example

# %%
class LittleMeta(type):
    def __new__(cls, clsname, superclasses, attributedict):
        print('clsname: ', clsname)
        print('superclasses: ', superclasses)
        print('attributedict: ', attributedict)
        return type.__new__(cls, clsname, superclasses, attributedict)


# %%
class S:
    pass

class A(S, metaclass=LittleMeta):
    pass


# %%
a = A()

# %%
print(a)

# %%
a.__class__.__bases__

# %%
a.__class__

# %%
type(S)


# %% [markdown]
# ## Factory Method #1

# %%

# %%
class ObjectFactory:
    def __init__(self):
        self._builders = {}

    def register_builder(self, key, builder):
        self._builders[key] = builder
        print(builder)

    def create(self, key, **kwargs):
        builder = self._builders.get(key)

        if not builder:
            raise ValueError(key)

        return builder(**kwargs)


# %%
class ServiceHYDAT:
    def __init__(self, db_hdl):
        self._db_hdl = db_hdl
        print(f'ServiceHYDAT.__init__(): dbl_hdl = {self._db_hdl}')

    def get_obs(self):
        print('Getting some data')
        pass


# %%
class ServiceHYDATBuilder:
    def __init__(self):
        print('ServiceHYDATBuilder.__init__()')
        self._instance = None

    def __call__(self, HYDAT_db_filename, **_ignored):
        print('ServiceHYDATBuilder.__call__()')
        if not self._instance:
            print('    *new* ServiceHYDAT instance')
            db_hdl = self.connect(HYDAT_db_filename)
            self._instance = ServiceHYDAT(db_hdl)
            print(self._instance)
        return self._instance

    def connect(self, db_filename):
        # This would connect to the sqlite3 db and return the handle
        print(f'Connecting to {db_filename}')
        return 'mydb_hdl'


# %%
factory = ObjectFactory()
factory.register_builder('HYDAT', ServiceHYDATBuilder())

# %%
config = {'HYDAT_db_filename': 'somefile.sqlite'}

hydat = factory.create('HYDAT', **config)

# %%
print(hydat)

# %%
hydat.get_obs()

# %%
id(hydat)

# %%
id(factory.create('HYDAT', **config))

# %%

# %% [markdown]
# ## Factory method approach #2

# %% language="javascript"
# IPython.notebook.kernel.restart()

# %%
from abc import ABCMeta, abstractmethod
from typing import Callable


# %%
class DataFactory:
    registry = {}

#     @classmethod
#     def register(cls, name: str) -> Callable:
#         def inner_wrapper(wrapped_class) -> Callable:
#             cls.registry[name] = wrapped_class
#             print(f'DataFactory: registered service {name}')
#             return wrapped_class
#         return inner_wrapper
    
    def register(self, name, connector):
        self.registry[name] = connector
        print(connector)
        
#     @classmethod
    def create_service(cls, name:str, **kwargs):
        svc_class = cls.registry.get(name)
        
        if not svc_class:
            raise ValueError(name)

        print(f'DataFactory.create_service(): {name} service retrieved')
        print(svc_class)
        return svc_class(**kwargs)



# %%
class ServiceBase(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def get_obs(self):
        # print('Getting some data')
        pass


# %%
# @DataFactory.register('HYDAT')
class ServiceHYDAT(ServiceBase):
    def __init__(self, db_hdl):
        self._db_hdl = db_hdl
        print(f'ServiceHYDAT.__init__(): dbl_hdl = {self._db_hdl}')
    
    def get_obs(self):
        print('Get HYDAT data')
        


# %%
# @DataFactory.register('HYDAT')
class ServiceHYDAT_connector():
    def __init__(self, **kwargs):
        print('ServiceHYDAT_connector.__init__()')
        self._instance = None

    def __call__(self, HYDAT_db_filename, **_ignored):
        print('ServiceHYDAT_connector.__call__()')
        if not self._instance:
            print('    *new* ServiceHYDAT instance')
            db_hdl = self.connect(HYDAT_db_filename)
            self._instance = ServiceHYDAT(db_hdl)
        return self._instance

    def connect(self, db_filename):
        # This would connect to the sqlite3 db and return the handle
        print(f'Connecting to {db_filename}')
        return 'mydb_hdl'    


# %%
dsources = DataFactory()
dsources.register('HYDAT', ServiceHYDAT_connector())

config = {'HYDAT_db_filename': 'somefile.sqlite'}
ds_hydat = dsources.create_service('HYDAT', **config)

# %%
ds_hydat

# %%
ds_hydat.get_obs()

# %%
DataFactory.registry

# %%
bb = ServiceHYDAT()

# %%
bb = ServiceHYDAT_connector('blah')

# %%
