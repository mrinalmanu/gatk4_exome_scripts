import sys
import os 
import fcntl
import time

def check_threads(batch_ID, task_ID,start_finish, threads_change, max_nr_threads):
	to_lock = open("pipelines/"+batch_ID+".CPUlog",'r+')    
	fcntl.flock(to_lock, fcntl.LOCK_EX)
	log=to_lock.read().splitlines()
	lastline= log[-1]
	
	threads_in_use=lastline.split()[-1]
	threads_update= int(threads_in_use) + int(threads_change)
	
	if threads_update > max_nr_threads:
		fcntl.flock(to_lock, fcntl.LOCK_UN)
		to_lock.close()
		return "wait"
		
	if int(threads_update) <= int(max_nr_threads):			
		time_stamp=time.time()
		to_lock.write(str(time_stamp)+" "+task_ID+" "+start_finish+" "+str(threads_update)+"\n")
		fcntl.flock(to_lock, fcntl.LOCK_UN)
		to_lock.close()
		return "go"
