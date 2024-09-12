# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 17:51:38 2024

@author: nesseler
"""

ePhys_parameter_dict = {'rheobase_abs' : {'ax_min' : -50,
                                          'ax_max' : 250,
                                          'pad' : 5,
                                          'step' : 50,
                                          'stepminor' : 10,
                                          'label' : 'Absolute rheobase [pA]'},
                        'rheobase_rel' : {'ax_min' : -50,
                                          'ax_max' : 250,
                                          'pad' : 5,
                                          'step' : 50,
                                          'stepminor' : 10,
                                          'label' : 'Relative rheobase [pA]'},
                        'rheobasespike_ttopeak' : {'ax_min' : 0.5,
                                                   'ax_max' : 1.5,
                                                   'pad' : (1.5 - 0.5) * 0.05,
                                                   'step' : 0.5,
                                                   'stepminor' : 0.1,
                                                   'label' : 'Rheobase spike time to peak [ms]'},
                        'rheobasespike_trise' : {'ax_min' : 0.2,
                                                 'ax_max' : 0.5,
                                                 'pad' : (0.5 - 0.2) * 0.05,
                                                 'step' : 0.1,
                                                 'stepminor' : 0.05,
                                                 'label' : 'Rheobase spike rise time [ms]'},
                        'rheobasespike_ttoAHP' : {'ax_min' : 1,
                                                  'ax_max' : 4.5,
                                                  'pad' : (4.5 - 1) * 0.05,
                                                  'step' : 1,
                                                  'stepminor' : 0.1,
                                                  'label' : 'Rheobase spike time to (fast) AHP [ms]'},
                        'rheobasespike_FWHM' : {'ax_min' : 0.5,
                                                'ax_max' : 2,
                                                'pad' : (2 - 0.5) * 0.05,
                                                'step' : 0.5,
                                                'stepminor' : 0.1,
                                                'label' : 'Rheobase spike FWHM [ms]'}}