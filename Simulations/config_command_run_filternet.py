#!/usr/bin/env python3
from bmtk.simulator import filternet
import nest
import sys

def main(config_file):
    """
    Put the config path/name argument in the python command
    """
    configure = filternet.Config.from_json(config_file[1])
    configure.build_env()

    network = filternet.FilterNetwork.from_config(configure)
    sim = filternet.FilterSimulator.from_config(configure, network)
    sim.run()

if __name__ == '__main__':
    main(sys.argv)
