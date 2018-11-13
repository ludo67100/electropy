#-*- coding: utf-8 -*-
"""
Created on Wed Oct 10 11:40:55 2018

@author: ludovic.spaeth
"""

'''
electroPy is a package using neo framework to analyze in vitro ephy recordings easily
'''

import logging

logging_handler = logging.StreamHandler()


from electroPy.Spike_Record import Spike_Record
from electroPy.WinWcpIo import WinWcpIo
from electroPy.SpikeDetection import SpikeDetection
from electroPy.HdF5IO import HdF5IO