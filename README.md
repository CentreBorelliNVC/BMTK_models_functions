To simulate l4_to_l4 network : 

1/Run ../Networks_builder/position_building_V1.py in order to build the V1 nodes and notably l4 nodes (nodes already saved in ../Networks/nodes if you wish to skip this step)

2/Run ../Networks_builder/l4_to_l4_connections.py in order to build the l4_to_l4 connections according to gathered data (edges already saved in ../Networks/l4_to_l4 if you whish to skip this step)

3/Run ../Simulations/run_pointnet.py in order to simulate the network according to ../Config_files/config.l4_10percent.json info (1000pA input on exc1 neurons)
