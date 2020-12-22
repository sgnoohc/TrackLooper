#!/bin/env python

import plottery_wrapper as p

p.dump_plot(
        fnames=[
            "debug.root",
            ],
        dirname="plots_multiplicity",
        extraoptions={
            "print_yield":True,
            }
        )
