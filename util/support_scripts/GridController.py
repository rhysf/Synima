# python LSF controllers
# based on Brian Haas' perl grid submission tools

import random, time
import os, sys, socket, subprocess, re

class GridController:
	"""
	To initialize an instance of the class:
	test_controller = BsubController.BsubController(test_cmd_list, ...)
	
	To run commands on the grid using that instance:
	test_controller.run_grid_submission()
	
	To get failed commands upon completion (retuns a list of tuples (command,job_id,return_value), no failed commands returns an empty list):
	test_controller.get_failed_cmds()
	
	To clean up logs:
	test_controller.clean_logs()
	
	"""
	def __init__(self, command_list, platform = 'LSF', queue = 'hour', dotkits = [], cmds_per_node = 50, memory = False, cpus=False, mount_test = False, max_nodes=500, debug = False, project = False):
		self.command_list = command_list
		self.queue = queue
		self.cmds_per_node = cmds_per_node
		self.memory = memory
		self.cpus = cpus
		self.mount_test = mount_test
		self.max_nodes = max_nodes
		self.nodes_in_progress = {}
		self.cmd_index_to_job_id = {}
		self.cmd_index_to_shell_script = {}
		self.job_id_to_submission_time = {}
		self.debug = debug  ## for debugging, enables extra messages
		self.project = project
		self.platform = platform
		self.dotkits = dotkits
		
		self.RESORT_TO_POLLING_TIME = 900
		
		self.num_cmds_launched = 0
		self.num_nodes_used = 0
		self.num_cmds = len(self.command_list)
		
		os.umask(0000)
		
		## make log directory
		
		self.log_id = str(random.randrange(10000000,100000000))
		self.log_dir_name = 'bsub.' + self.log_id
		if os.path.exists(self.log_dir_name):
			log_except_str = 'log_dir ' + self.log_dir_name + ' already exists'
			raise Exception(log_except_str)
		else:
			os.makedirs(self.log_dir_name)
			
		## write commands to log directory	
		
		file = open((self.log_dir_name + '/cmds_list.txt'), 'w')
		command_index = 0
		for command in self.command_list:
			file.write('index(' + str(command_index) + ')\t' + command + '\n')
			command_index += 1
		file.close()
		
		## finish logging setup
		
		self.command_dir = self.log_dir_name + '/cmds'
		self.retvals_dir = self.log_dir_name + '/retvals'
		self.monitor_dir = self.log_dir_name + '/monitor'
		for dir in [self.command_dir, self.retvals_dir, self.monitor_dir]:
			os.makedirs(dir)
		
	def get_command_list(self):
		return self.command_list
	
	def get_log_dir_name(self):
		return self.log_dir_name
		
	def write_pid_file(self):
		hostname = socket.gethostname()		
		file = open((self.log_dir_name + '/' + hostname + '.pid'), 'w')
		file.write(str(os.getpid()))
		file.close()
		
	def run_grid_submission(self):
		self.write_pid_file()

		while self.num_cmds_launched < self.num_cmds:
 			self.num_cmds_launched = self.submit_job()
 			self.num_nodes_used = self.get_num_nodes_in_use()
 			sys.stdout.write('\rCMDS: ' + str(self.num_cmds_launched) + '/' + str(self.num_cmds) + '  [' + str(self.num_nodes_used) + '/' + str(self.max_nodes) + ' nodes in use]     ')
			sys.stdout.flush()
			if self.num_nodes_used >= self.max_nodes:
				num_nodes_finished = self.wait_for_completions()
				self.num_nodes_used = self.num_nodes_used - num_nodes_finished
		
		print('\nAll cmds submitted to grid.  Now waiting for them to finish.')
		## wait for rest to finish
		num_nodes_finished = self.wait_for_completions()
		while num_nodes_finished:
			num_nodes_finished = self.wait_for_completions()
			self.num_nodes_used = self.num_nodes_used - num_nodes_finished
			sys.stdout.write('\rCMDS: ' + str(self.num_cmds_launched) + '/' + str(self.num_cmds) + '  [' + str(self.num_nodes_used) + '/' + str(self.max_nodes) + ' nodes in use]     ')
			sys.stdout.flush()
		print('\nAll nodes completed.  Now auditing job completion status values.')

		self.get_exit_values()
		
		num_successes = 0
		num_failures = 0
		num_unknown = 0
		
		for retval in self.retvals:
			try:
				int(retval)
				if int(retval) == 0:
					num_successes = num_successes + 1
				else:
					num_failures = num_failures + 1
			except:
				num_unknown = num_unknown + 1

		self.write_result_summary(num_successes, num_failures, num_unknown)
		
		if num_successes == self.num_cmds:
			print('All commands completed successfully.')
		else:
			print('num_success: ' + str(num_successes) + ' num_fail: ' + str(num_failures) + ' num_unknown: ' + str(num_unknown))

		print('Finished')   	
		
	def submit_job(self):
		orig_num_cmds_launched = self.num_cmds_launched
		
		shell_script = self.command_dir + '/S' + str(self.log_id) + '.' + str(self.num_cmds_launched) + '.sh'
		file = open(shell_script, 'w')	
		file.write('#!/bin/sh\n\n')
		file.write('## add any special environment settings\n\n')
		file.write('echo HOST: $HOSTNAME\n')
		file.write('echo HOST: $HOSTNAME >&2\n\n')
		if self.platform == 'GridEngine' or self.platform == 'UGER' or self.dotkits:
			file.write('source /broad/software/scripts/useuse\n')
		if self.platform == 'GridEngine':
			file.write('reuse GridEngine8\n')
		if self.platform == 'UGER':
			file.write('reuse UGER\n')		
		if self.dotkits:
			for dotkit in self.dotkits:
				file.write('reuse ' + dotkit + '\n')
		
		num_cmds_written = 0
		monitor_started = self.monitor_dir + '/' + str(self.num_cmds_launched) + '.started'
		monitor_finished = self.monitor_dir + '/' + str(self.num_cmds_launched) + '.finished'
		
		cmd_indices_prepped = []
		
		while (self.num_cmds_launched < self.num_cmds) and (num_cmds_written < self.cmds_per_node):  ## Brian's code had && instead of and
			next_cmd_index = self.num_cmds_launched
			cmd_string = self.command_list[next_cmd_index]
			self.cmd_index_to_shell_script[next_cmd_index] = shell_script
			cmd_indices_prepped.append(next_cmd_index)
			
			retval_bin = str(next_cmd_index / 1000)
			retval_subdir = self.retvals_dir + '/' + retval_bin
			if not os.path.exists(retval_subdir):
				os.makedirs(retval_subdir)
			
			file.write('## Command index ' + str(next_cmd_index) + '\n')
			file.write('touch ' + monitor_started + '\n')
			file.write(cmd_string + '\n')
			file.write('echo $? >> ' + retval_subdir + '/entry_' + str(next_cmd_index) + '.ret\n\n')
			
			self.num_cmds_launched += 1
			num_cmds_written += 1

		file.write('\nrm -f ' + monitor_started + '\n')
		file.write('touch ' + monitor_finished + '\n\n')
		file.write('exit 0\n\n')
		file.close()

		os.chmod(shell_script, 0775)
		
		if self.debug:
			print('Submitting ' + shell_script + ' to grid')
			
		script_basename = os.path.basename(shell_script)
		cmd = ''
		if self.platform == 'LSF':
			cmd = 'bsub -q ' + self.queue + ' -e ' + shell_script + '.stderr -o ' + shell_script + '.stdout'

			if self.memory:
				cmd = cmd + ' -R \"rusage[mem=' + str(self.memory) + ']\"'
				
			if self.cpus:
				cmd = cmd + ' -n ' + str(self.cpus) + ' -R \"span[hosts=1]\"'
		
			if self.queue == 'hour':
				cmd = cmd + ' -W 4:0'
				
			if self.project:
				cmd = cmd + ' -P ' + self.project
			
			if self.mount_test:
				cmd = cmd + ' -E \"/broad/tools/NoArch/pkgs/local/checkmount ' + self.mount_test + ' && [ -e ' + self.mount_test + ' ]\"'

		
		elif self.platform == 'GridEngine' or self.platform == 'UGER':
			cmd = 'qsub -V -cwd -q ' + self.queue + ' -e ' + shell_script + '.stderr -o ' + shell_script + '.stdout'

			if self.platform == 'GridEngine':
				if self.memory:
					cmd = cmd + ' -l h_vmem=' + str(self.memory) + 'G'
				
				if self.cpus:
					cmd = cmd + ' -pe smp_pe ' + str(self.cpus)
					
			elif self.platform == 'UGER':
				if self.memory:
					memory_setting = self.memory
					if self.cpus:
						memory_setting = int(self.memory/self.cpus) + (self.memory%self.cpus > 0)
					cmd = cmd + ' -l m_mem_free=' + str(memory_setting) + 'g'
				
				if self.cpus:
					cmd = cmd + ' -pe smp ' + str(self.cpus)
									
			if self.project:
				cmd = cmd + ' -P ' + self.project
			
			if self.mount_test:
				print('Mount test unavailable through GridEngine/UGER')
						
		else:
			raise Exception(('Unsupported grid platform: ' + self.platform))

		cmd = cmd + ' ' + shell_script + ' 2>&1'
		
		if self.debug:
			print(cmd)
 
 		#submission_return = subprocess.call(cmd, shell=True)
 		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
 		submission_out = process.communicate()[0]
 		submission_return = process.returncode
 		
 		if submission_return:
 			print('Grid failed to accept job: ' + cmd + '\n (ret ' + str(submission_return) + ')\n')
 			os.unlink(shell_script) # cleanup, try again later
 			time.sleep(120) # sleep 2 minutes for now.  Give the system time to recuperate if a problem exists
 			return orig_num_cmds_launched
 		else:
 			shell_script = os.path.basename(shell_script)
 			file = open((self.log_dir_name + '/job_ids.txt'), 'a')
 			job_pattern = re.compile(r'Job \<(\d+)\>')
 			if self.platform == 'GridEngine' or self.platform == 'UGER':
 				job_pattern = re.compile(r'Your job (\d+)')
 			matched = job_pattern.search(submission_out)
 			if matched:
 				job_id = matched.group(1)
 				file.write(job_id + '\t' + shell_script + '\n')
 				self.nodes_in_progress[monitor_finished] = job_id ## hope this is right
 				for cmd_index in cmd_indices_prepped:
 					self.cmd_index_to_job_id[cmd_index] = job_id
 				if self.debug:
 					print('job id: ' + str(job_id) + ' cmd index: ' + str(cmd_index))
 				self.job_id_to_submission_time[job_id] = int(time.time())
 				## self.job_id_tester = job_id ## for testing only
			else:
				raise Exception(('Fatal error, couldn\'t extract Job ID from submission text: ' + job_id_text))
			file.close()
			#time.sleep(15) # wait just a short while to give the system a few seconds to act on the submitted jobs
			return(self.num_cmds_launched)

	
	def get_num_nodes_in_use(self):
		num_nodes_in_use = len(self.nodes_in_progress.keys())
		if self.debug:
			print('Num nodes currently in use: ' + str(num_nodes_in_use))
		return num_nodes_in_use
		
	def wait_for_completions(self):
		
		if self.debug:
			print('Running wait_for_completions()')
		
		seen_finished = 0
		done = []
		
		while not seen_finished:
			## check to see if there are any jobs remaining:
			if self.get_num_nodes_in_use() == 0:
				## no jobs in the queue
				if self.debug:
					print('No nodes in use, exiting wait')
				return 0
				
			## check for finished jobs
			for monitor_file in self.nodes_in_progress.keys():
				if os.path.isfile(monitor_file):
					done.append(monitor_file)
					seen_finished = 1
				else:
					job_id = self.nodes_in_progress[monitor_file]
					time_launched = self.job_id_to_submission_time[job_id]
					current_time = int(time.time())
					## see if an hour has passed
					if (current_time - time_launched) >= self.RESORT_TO_POLLING_TIME:
						## poll the system directly:
						if not self.is_job_running_or_pending_on_grid(job_id):
							done.append(monitor_file)
							seen_finished = 1
						else:
							## reset submission time to delay next polling time
							self.job_id_to_submission_time[job_id] = int(time.time())
					
			if seen_finished:
				for monitor_file in done:
					job_id = self.nodes_in_progress[monitor_file]
					if self.debug:
						print('job[' + str(job_id) + ']: ' + str(monitor_file) + ' is finished')
					del self.nodes_in_progress[monitor_file]
					del self.job_id_to_submission_time[job_id]
				return len(done)
			else:
				## wait a while and check again
				if self.debug:
					print('Waiting for jobs to finish')
				time.sleep(15)

	
	def is_job_running_or_pending_on_grid(self, job_id):
		## job_id = self.job_id_tester ## for testing only
		if (int(time.time()) - self.job_id_to_submission_time[job_id]) < self.RESORT_TO_POLLING_TIME:
			return('TOO_SOON')
		
		if self.debug:
			print('Polling grid to check status of job: ' + str(job_id))

		attempts = 0		
		while attempts < 5:
			cmd = 'bjobs ' + str(job_id)
 			if self.platform == 'GridEngine' or self.platform == 'UGER':
 				cmd = 'qstat -s za'
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
			submission_out = process.communicate()[0]
			submission_return = process.returncode
		
			if not submission_return:
				if self.debug:
					print('Submission out: ' + submission_out)
				split_out = submission_out.split('\n')
				for split_out_line in split_out:
					split_out2 = split_out_line.split()
					try:
						test_index = split_out2[0]
						if split_out2[0] == job_id:
							state = split_out2[2]
							if self.platform == 'GridEngine' or self.platform == 'UGER':
								state = split_out2[4]
							if state == 'DONE' or state == 'EXIT' or state == 'UNKWN' or state == 'ZOMBI' or state == 'z' or state == 'Eqw':
								return 0
							else:
								self.job_id_to_submission_time[job_id] = int(time.time())
								return state
					except:
						return 0 ### job missing from qstat output
							
			attempts = attempts + 1
			time.sleep(15)
 					
		print('No record of job_id ' + str(job_id) + ', setting as state unknown\n')
		return 0 ## Brian returns that as unknown, but it results in the same as UNKWN
				
    
		
	def get_ret_filename(self, cmd_index):
		retval_bin = str(cmd_index/1000);
		retval_file = self.retvals_dir + '/' + retval_bin + '/entry_' + str(cmd_index) + '.ret'
		return retval_file
		
	def clean_logs(self):
		pass
		cmd = 'rm -rf ' + self.log_dir_name
		return_code = subprocess.call(cmd, shell=True)
		return return_code
		
	def write_result_summary(self, num_successes, num_failures, num_unknown):
		status = 'failure'
		if num_failures == 0 and num_unknown == 0:
			status = 'success'
		
		file = open((self.log_dir_name + '/bsub.finished.' + status), 'a')
		file.write('num_successes: ' + str(num_successes) + '\n')
		file.write('num_failures: ' + str(num_failures) + '\n')
		file.write('num_unknown: ' + str(num_unknown) + '\n')
		file.close()
		
	def get_failed_cmds(self):
		failed_cmds = []
		for i in range(len(self.retvals)):
			if self.retvals[i]:
				failed_cmds.append((self.command_list[i],self.cmd_index_to_job_id[i],self.retvals[i],self.cmd_index_to_shell_script[i]))
		return failed_cmds
		
	def get_exit_values(self):
		self.retvals = []
		if self.debug:
			print('Processing ' + self.retvals_dir)
		for i in range(self.num_cmds):
			retval_file = self.get_ret_filename(i)
				
			if self.debug:
				print('file: ' + retval_file)
			try:
				os.path.getsize(retval_file)
				file = open(retval_file, 'r')
				retval_string = file.read()
				retval_string = ''.join(retval_string.split())
				if self.debug:
					print('retval: ' + retval_string)
				self.retvals.append(int(retval_string))
				file.close()
			except:
				self.retvals.append('FILE_NOT_EXISTS')


