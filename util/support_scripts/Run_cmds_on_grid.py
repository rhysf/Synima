#!/usr/bin/env python
# LSF/GridEngine/UGER grid running wrapper
# based on Brian Haas' perl grid submission tools

import GridController
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('command_list', help='list of commands to run', type=str)
parser.add_argument('--platform', help='LSF|GridEngine|UGER (default = UGER)', type=str)
parser.add_argument('--queue', help='queue to submit too (default = short)', type=str)
parser.add_argument('--mem', help='memory in G (default is unspecified)', type=int)
parser.add_argument('--cpus', help='cpus (default is unspecified)', type=int)
parser.add_argument('--cmds_per_node', help='commands for each grid node to process (default=50)', type=int)
parser.add_argument('--throttle_nodes', help='max number of commands to run simultaneously (default=500)', type=int)
parser.add_argument('--mount_test', help='directory for grid nodes to verify proper mounting (LSF only)', action='store_true')
parser.add_argument('--dotkits', help='list of dotkits to load', type=str)

args = parser.parse_args()

platform = 'UGER'
queue = 'short'
mem = False
cpus = False
cmds_per_node = 50
max_nodes = 500
mount_test = args.mount_test
dotkits = False

if args.platform:
	platform = args.platform
if args.queue:
	queue = args.queue
if args.mem:
	mem = args.mem
if args.cpus:
	cpus = args.cpus
if args.cmds_per_node:
	cmds_per_node = args.cmds_per_node
if args.throttle_nodes:
	max_nodes = args.throttle_nodes
if args.dotkits:
	dotkits = args.dotkits
	
if platform != 'LSF' and platform != 'GridEngine' and platform != 'UGER':
    raise Exception("Error, LSF, GridEngine, and UGER are the only recognized platforms.\n")	

if (platform == 'GridEngine' or platform == 'UGER')	and mount_test:
	print("Warning, mount test option is not passed through GridEngine/UGER.\n")


cmds_list = []
file = open(args.command_list, 'r')
for line in file:
		cmds_list.append(line.strip())
file.close()

dotkits_list = []
if dotkits:
	file = open(dotkits, 'r')
	for line in file:
		dotkits_list.append(line.strip())
	file.close()

controller = GridController.GridController(cmds_list, platform=platform, queue=queue, memory=mem, cpus=cpus, cmds_per_node=cmds_per_node, max_nodes=max_nodes, mount_test=mount_test, dotkits=dotkits_list)

controller.run_grid_submission()

failed_cmds = controller.get_failed_cmds()
if failed_cmds:
	for failed in failed_cmds:
		print(failed[0] + ' in job ' + str(failed[1]) + ' failed with ret ' + str(failed[2]) + ', see ' + str(failed[3]) + '.stdout and ' + str(failed[3]) + '.stderr')
else:
	print('No failed commands.')
	controller.clean_logs()



