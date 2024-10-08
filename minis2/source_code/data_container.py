import numpy as np
import pyabf

class DataContainer:
  """
  Minis data container encapsulating voltage and current data traces.
  
  Attributes:
    voltage (ndarray): a shape-(1, n) ndarray with a membrane potential trace
      recording.
    current (ndarray): a shape-(1, n) ndarray with an injected transmembrane
      current trace recording.
    filename (str): a shape-(1, m) data filename string.
    state (str): a shape-(1, l) data derivation state string: 'passed in' or
      'loaded'.
  """

  def __init__(self,
               voltage: np.ndarray = np.empty((0,0)),
               current: np.ndarray = np.empty((0,0)),
               filename: str = ''):
    """
    Initialises a data DataContainer instance.
    
    Args:
      voltage (ndarray, optional): a shape-(1, n) ndarray with a membrane
        potential trace recording.
      current (ndarray, optional): a shape-(1, n) ndarray with an injected
        transmembrane current trace recording.
      filename (str, optional): a shape-(1, m) data filename string. If
        supplied, supplied voltage and current data are overriden by the data
        loaded from the file. 
    """

    self.voltage = voltage
    self.current = current
    self.filename = filename
    self.state = 'passed in'

    # Load data from the file, if a nonempty filename is supplied
    if len(filename):
      abf = pyabf.ABF(filename)
      self.voltage = abf.sweepY
      self.current = abf.sweepC
      self.filename = abf.abfFilePath
      self.state = 'loaded'


minis_data = DataContainer(filename='p103a_0100-0104.abf')
print(minis_data.voltage)
print(minis_data.current)
print(minis_data.filename)
print(minis_data.state)