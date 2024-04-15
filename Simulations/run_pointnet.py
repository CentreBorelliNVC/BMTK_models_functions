#!/usr/bin/env python3
from bmtk.simulator import pointnet
import nest
from bmtk.analyzer.spike_trains import plot_raster
from bmtk.analyzer.compartment import plot_traces

def main(config_file):
    configure = pointnet.Config.from_json(config_file)
    configure.build_env()

    network = pointnet.PointNetwork.from_config(configure)
    sim = pointnet.PointSimulator.from_config(configure, network)
    sim.run()
    #pointnet.nrn.quit_execution()
    _ = plot_raster(config_file=config_file, group_by='pop_name', show=True)
    _ = plot_traces(config_file=config_file,report_name='V_m_report')


if __name__ == '__main__':
    main('../Config_files/config.l4_10percent.json')
