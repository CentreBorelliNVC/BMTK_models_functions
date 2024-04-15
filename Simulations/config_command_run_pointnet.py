#!/usr/bin/env python3
from bmtk.simulator import pointnet
import nest
from bmtk.analyzer.spike_trains import plot_raster
from bmtk.analyzer.compartment import plot_traces
import sys

def main(config_file):
    """
    Put the config path/name argument in the python command
    """
    configure = pointnet.Config.from_json(config_file[1])
    configure.build_env()

    network = pointnet.PointNetwork.from_config(configure)
    sim = pointnet.PointSimulator.from_config(configure, network)
    sim.run()
    #pointnet.nrn.quit_execution()
    _ = plot_raster(config_file=config_file, group_by='pop_name', show=True)
    _ = plot_traces(config_file=config_file,report_name='V_m_report')


if __name__ == '__main__':
    main(sys.argv)
