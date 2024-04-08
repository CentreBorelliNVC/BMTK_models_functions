from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder


factor=0.1

#L1  : 
filename='../Components/nrrd_files/structure_593.nrrd'
l1=CellLocations('l1')
l1.dmin=14.8
l1.CCF_orientation=True
max_dens_all=43199.7*factor
l1.add_positions_nrrd(filename,max_dens_all,pop_names=['htr3a1','htr3a2'],partitions=[0.34,0.66],method='prog',verbose=True)

#print('htr3a : ',len(l1.htr3a1.positions),len(l1.htr3a2.positions))


net = NetworkBuilder('cortex_l1_'+str(factor))
net.add_nodes(N=l1.htr3a1.N,
              positions=l1.htr3a1.positions,
              pop_name='htr3a1', location='VisL1', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L1_Htr3a1_529883902_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l1.htr3a2.N,
              positions=l1.htr3a2.positions,
              pop_name='htr3a2', location='VisL1', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L1_htr3a2_566307069_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )             
net.build()
net.save_nodes(output_dir='../Networks/nodes')


#L2/3  : 
filename='../Components/nrrd_files/structure_821.nrrd'
l23=CellLocations('l23')
l23.dmin=13.8
l23.CCF_orientation=True
max_dens_all=108787.9*factor
l23.add_positions_nrrd(filename,max_dens_all,pop_names=['exc','PV1','PV2','PV3','PV4','SST1','SST2','SST3','VIP1','VIP2','VIP3','VIP4','VIP5','htr3a1','htr3a2'],partitions=[0.75,0.01,0.01,0.01,0.01,0.00667,0.00667,0.00666,0.01,0.01,0.01,0.01,0.01,0.07,0.07],method='prog',verbose=True)

#print(len(l23.exc.positions),len(l23.PV1.positions),len(l23.PV2.positions),len(l23.PV3.positions),len(l23.PV4.positions),len(l23.SST1.positions),len(l23.SST2.positions),len(l23.SST3.positions),len(l23.VIP1.positions),len(l23.VIP2.positions),len(l23.VIP3.positions),len(l23.VIP4.positions),len(l23.VIP5.positions),len(l23.htr3a1.positions),len(l23.htr3a2.positions))


net = NetworkBuilder('cortex_l23_'+str(factor))
net.add_nodes(N=l23.exc.N,
              positions=l23.exc.positions,
              pop_name='exc', location='Visl23', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_exc_515020572_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.PV1.N,
              positions=l23.PV1.positions,
              pop_name='PV1', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_pvalb1_566290552_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l23.PV2.N,
              positions=l23.PV2.positions,
              pop_name='PV2', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_pvalb2_489931670_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.PV3.N,
              positions=l23.PV3.positions,
              pop_name='PV3', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_pvalb3_587846526_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l23.PV4.N,
              positions=l23.PV4.positions,
              pop_name='PV4', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_pvalb4_591765645_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.SST1.N,
              positions=l23.SST1.positions,
              pop_name='SST1', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_sst1_566284338_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.SST2.N,
              positions=l23.SST2.positions,
              pop_name='SST2', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_sst2_561761257_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.SST3.N,
              positions=l23.SST3.positions,
              pop_name='SST3', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_sst3_566316602_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.VIP1.N,
              positions=l23.VIP1.positions,
              pop_name='VIP1', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_vip1_571019402_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l23.VIP2.N,
              positions=l23.VIP2.positions,
              pop_name='VIP2', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_vip2_566285504_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.VIP3.N,
              positions=l23.VIP3.positions,
              pop_name='VIP3', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_vip3_566729508_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l23.VIP4.N,
              positions=l23.VIP4.positions,
              pop_name='VIP4', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_vip4_587848209_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.VIP5.N,
              positions=l23.VIP5.positions,
              pop_name='VIP5', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_vip5_591774123_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )          
net.add_nodes(N=l23.htr3a1.N,
              positions=l23.htr3a1.positions,
              pop_name='htr3a1', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_htr3a1_591248063_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l23.htr3a2.N,
              positions=l23.htr3a2.positions,
              pop_name='htr3a2', location='Visl23', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L23_htr3a2_557414710_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )             
net.build()
net.save_nodes(output_dir='../Networks/nodes')

#L4 :
filename='../Components/nrrd_files/structure_721.nrrd'
l4=CellLocations('l4')
l4.dmin=13.2
l4.CCF_orientation=True
max_dens_all=130640.6*factor
l4.add_positions_nrrd(filename,max_dens_all,pop_names=['exc1','exc2','exc3','exc4','exc5','exc6','exc7','PV1','PV2','SST1','SST2','SST3','VIP1','VIP2','VIP3','VIP4','htr3a'],partitions=[0.107,0.107,0.107,0.107,0.107,0.107,0.108,0.045,0.045,0.01667,0.01667,0.01666,0.005,0.005,0.005,0.005,0.09],method='prog',verbose=True)

#print(len(l4.exc1.positions),len(l4.exc2.positions),len(l4.exc3.positions),len(l4.exc4.positions),len(l4.exc5.positions),len(l4.exc6.positions),len(l4.exc7.positions),len(l4.PV1.positions),len(l4.PV2.positions),len(l4.SST1.positions),len(l4.SST2.positions),len(l4.SST3.positions),len(l4.VIP1.positions),len(l4.VIP2.positions),len(l4.VIP3.positions),len(l4.VIP4.positions),len(l4.htr3a.positions))


net = NetworkBuilder('cortex_l4_'+str(factor))
net.add_nodes(N=l4.exc1.N,
              positions=l4.exc1.positions,
              pop_name='exc1', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc1_566290754_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.exc2.N,
              positions=l4.exc2.positions,
              pop_name='exc2', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc2_566374978_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.exc3.N,
              positions=l4.exc3.positions,
              pop_name='exc3', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc3_512192761_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.exc4.N,
              positions=l4.exc4.positions,
              pop_name='exc4', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc4_566304492_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.exc5.N,
              positions=l4.exc5.positions,
              pop_name='exc5', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc5_566293986_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.exc6.N,
              positions=l4.exc6.positions,
              pop_name='exc6', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc6_517842022_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.exc7.N,
              positions=l4.exc7.positions,
              pop_name='exc7', location='Visl4', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_exc7_566290186_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.PV1.N,
              positions=l4.PV1.positions,
              pop_name='PV1', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_pvalb1_485694390_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l4.PV2.N,
              positions=l4.PV2.positions,
              pop_name='PV2', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_pvalb2_489931670_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.SST1.N,
              positions=l4.SST1.positions,
              pop_name='SST1', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_sst1_566729393_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.SST2.N,
              positions=l4.SST2.positions,
              pop_name='SST2', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_sst2_566737689_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.SST3.N,
              positions=l4.SST3.positions,
              pop_name='SST3', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_sst3_587845733_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.VIP1.N,
              positions=l4.VIP1.positions,
              pop_name='VIP1', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_vip1_566305624_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l4.VIP2.N,
              positions=l4.VIP2.positions,
              pop_name='VIP2', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_vip2_562006495_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l4.VIP3.N,
              positions=l4.VIP3.positions,
              pop_name='VIP3', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_vip3_587848241_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l4.VIP4.N,
              positions=l4.VIP4.positions,
              pop_name='VIP4', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_vip4_591762171_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )         
net.add_nodes(N=l4.htr3a.N,
              positions=l4.htr3a.positions,
              pop_name='htr3a', location='Visl4', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L4_ht3ar_561769814_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )           
net.build()
net.save_nodes(output_dir='../Networks/nodes')

#L5 : 
filename='../Components/nrrd_files/structure_778.nrrd'
l5=CellLocations('l5')
l5.dmin=13.4
l5.CCF_orientation=True
max_dens_all=105268.4*factor
l5.add_positions_nrrd(filename,max_dens_all,pop_names=['CF','IT1','IT2','PV','SST1','SST2','SST3','SST4','VIP','htr3a'],partitions=[0.24,0.24,0.24,0.12,0.02,0.02,0.02,0.02,0.02,0.06],method='prog',verbose=True)

#print(len(l5.CF.positions),len(l5.IT1.positions),len(l5.IT2.positions),len(l5.PV.positions),len(l5.SST1.positions),len(l5.SST2.positions),len(l5.SST3.positions),len(l5.SST4.positions),len(l5.VIP.positions),len(l5.htr3a.positions))


net = NetworkBuilder('cortex_l5_'+str(factor))
net.add_nodes(N=l5.CF.N,
              positions=l5.CF.positions,
              pop_name='CF', location='Visl5', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_CF_486556800_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l5.IT1.N,
              positions=l5.IT1.positions,
              pop_name='IT1', location='Visl5', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_IT1_557185822_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l5.IT2.N,
              positions=l5.IT2.positions,
              pop_name='IT2', location='Visl5', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_IT2_486557270_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l5.PV.N,
              positions=l5.PV.positions,
              pop_name='PV', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_pvalb_512521589_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l5.SST1.N,
              positions=l5.SST1.positions,
              pop_name='SST1', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_sst1_520598286_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l5.SST2.N,
              positions=l5.SST2.positions,
              pop_name='SST2', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_sst2_550252140_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l5.SST3.N,
              positions=l5.SST3.positions,
              pop_name='SST3', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_sst3_569476228_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l5.SST4.N,
              positions=l5.SST4.positions,
              pop_name='SST4', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_sst4_587864555_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )

net.add_nodes(N=l5.VIP.N,
              positions=l5.VIP.positions,
              pop_name='VIP', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_vip_562006495_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )        
net.add_nodes(N=l5.htr3a.N,
              positions=l5.htr3a.positions,
              pop_name='htr3a', location='Visl5', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L5_htr3a_561769810_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )           
net.build()
net.save_nodes(output_dir='../Networks/nodes')

#L6 : 
filename='../Components/nrrd_files/structure_33.nrrd'
l6=CellLocations('l6')
l6.dmin=12.6
l6.CCF_orientation=True
max_dens_all=102263.2*factor
l6.add_positions_nrrd(filename,max_dens_all,pop_names=['exc1','exc2','exc3','PV1','PV2','PV3','PV4','SST1','SST2','SST3','VIP','htr3a'],partitions=[0.27,0.27,0.27,0.0175,0.0175,0.0175,0.0175,0.01667,0.01667,0.01666,0.01,0.06],method='prog',verbose=True)

#print(len(l6.exc1.positions),len(l6.exc2.positions),len(l6.exc3.positions),len(l6.PV1.positions),len(l6.PV2.positions),len(l6.PV3.positions),len(l6.PV4.positions),len(l6.SST1.positions),len(l6.SST2.positions),len(l6.SST3.positions),len(l6.VIP.positions),len(l6.htr3a.positions))


net = NetworkBuilder('cortex_l6_'+str(factor))
net.add_nodes(N=l6.exc1.N,
              positions=l6.exc1.positions,
              pop_name='exc1', location='Visl6', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_exc1_566303013_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.exc2.N,
              positions=l6.exc2.positions,
              pop_name='exc2', location='Visl6', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_exc2_587864254_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.exc3.N,
              positions=l6.exc3.positions,
              pop_name='exc3', location='Visl6', ei='e',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_exc3_565418574_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.PV1.N,
              positions=l6.PV1.positions,
              pop_name='PV1', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_pvalb1_587845769_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.PV2.N,
              positions=l6.PV2.positions,
              pop_name='PV2', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_pvalb2_573430548_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.PV3.N,
              positions=l6.PV3.positions,
              pop_name='PV3', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_pvalb3_566357258_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.PV4.N,
              positions=l6.PV4.positions,
              pop_name='PV4', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_pvalb4_592920664_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.SST1.N,
              positions=l6.SST1.positions,
              pop_name='SST1', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_sst1_587845932_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.SST2.N,
              positions=l6.SST2.positions,
              pop_name='SST2', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_sst2_538758825_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.SST3.N,
              positions=l6.SST3.positions,
              pop_name='SST3', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_sst3_520598286_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )
net.add_nodes(N=l6.VIP.N,
              positions=l6.VIP.positions,
              pop_name='VIP', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_vip_566736311_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )        
net.add_nodes(N=l6.htr3a.N,
              positions=l6.htr3a.positions,
              pop_name='htr3a', location='Visl6', ei='i',
              model_type='point_process',
              model_template='nest:glif_lif_asc_psc',
              dynamics_params='L6_htr3a_591407725_neuron_config.json'  # File containing glif_lif_asc_psc model parameters
             )           
net.build()
net.save_nodes(output_dir='../Networks/nodes')












