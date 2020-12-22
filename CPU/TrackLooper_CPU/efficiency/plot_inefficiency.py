#!/bin/env python

import plottery_wrapper as p

p.dump_plot(
        fnames=[
            "inefficiency.root",
            ],
        dirname="plots_inefficiency",
        extraoptions={
            "print_yield":True,
            }
        )
