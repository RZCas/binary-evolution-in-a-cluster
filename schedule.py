from mpi4py import MPI
import numpy as np
import time
from binary_evolution_with_flybys import inputParameters, evolve_binary, approximation_test, a_h
from scipy.stats import loguniform
# folder = '/nfs/st01/hpc-astro-gio10/ar2094/a_dependence_3/'
folder = 'output/a_dependence_3/'

def genTasks(numTasks):

    # np.random.seed(1000)  # run the same set of timed tasks

    t = 1.4e10
    output_file = ''

    # Inner binary parameters
    ecc = 0.5              # Eccentricity
    inc = 1           # Inclination with respect to the z-axis
    long_asc = 0            # Longitude of the ascending node
    arg_peri = 1.5    # Arugment of pericentre
    m1 = 5
    m2 = 5

    # Outer binary parameters
    ecc_out = 0.0         # Outer orbit eccentricity
    inc_out = 0.5             # Outer orbit inclination
    a_out = 0.5        # Outer semi-major axis in pc

    a_in_min = 0.1*a_h(m1,m2,a_out)              # Semi-major axis in AU
    a_in_max = 10*a_h(m1,m2,a_out)
    tmax = 3000000*60*60

    tasks=[]
    for i in range (numTasks):
        # a_in = (a_in_max-a_in_min)*np.random.random_sample()+a_in_min
        # a_in = loguniform.rvs (a_in_min, a_in_max)
        a_in = a_in_min*(a_in_max/a_in_min)**((i)/numTasks)
        #a_in = 1
        tasks.append(inputParameters(t=t, a_out=a_out, e_out=ecc_out, inc_out=inc_out, m1=m1, m2=m2, a=a_in, e=ecc, i=inc, Omega=long_asc, omega=arg_peri, output_file=output_file, forcePrecise=False, tmax=tmax, resume=True, includeWeakEncounters=True))

    return np.array(tasks)

# tags that can be applied to messages
WORKTAG = 1
STOPTAG = 2

def my_program():
    comm = MPI.COMM_WORLD
    id = comm.Get_rank()            #number of the process
    numProc = comm.Get_size()  #total number of processes

    if (id == 0) :
        numTasks = 20#(numProc-1)*1 # avg tasks per worker process
        inputs = genTasks(numTasks)
        #print(numTasks)
        #print(inputs, flush=True)
        superwise_work(inputs, comm, numProc)
    else:
        do_work(comm)

def superwise_work(inputs, comm, numProc):
    totalWork = inputs.size
    workcount = 0
    recvcount = 0
    # send out the first tasks to all workers
    for id in range(1, numProc):
        if workcount < totalWork:
            work=inputs[workcount]
            work.output_file = folder+str(workcount)+'.txt'
            work.output_file_2 = folder+'debug'+str(workcount)+'.txt'    
            comm.send(work, dest=id, tag=WORKTAG)
            workcount += 1
            print("master sent {} to {}".format(work, id), flush=True)

    # while there is still work,
    # receive result/requests for more work from a worker
    while (workcount < totalWork) :
        # receive next finished result
        stat = MPI.Status()
        result = comm.recv(source=MPI.ANY_SOURCE, status=stat)
        recvcount += 1
        workerId = stat.Get_source()
        print("master received {} from {}".format(result, workerId), flush=True)
        #send next work
        work=inputs[workcount]
        work.output_file = folder+str(workcount)+'.txt'
        work.output_file_2 = folder+'debug'+str(workcount)+'.txt'       
        comm.send(work, dest=workerId, tag=WORKTAG)
        workcount += 1
        print("master sent {} to {}".format(work, workerId), flush=True)

    # Receive results for outstanding work requests.
    while (recvcount < totalWork):
        stat = MPI.Status()
        result = comm.recv(source=MPI.ANY_SOURCE, status=stat)
        recvcount += 1
        workerId = stat.Get_source()
        print("master received {} from {}".format(result, workerId), flush=True)

    # Tell all workers to stop
    for id in range(1, numProc):
        comm.send(inputs[0], dest=id, tag=STOPTAG)


def do_work(comm):
    # loop till the tag says to stop
    while(True):
        stat = MPI.Status()
        input = comm.recv(source=0, tag=MPI.ANY_TAG, status=stat)
        # print("worker {} got {}".format(comm.Get_rank(), waitTime), flush=True)
        if (stat.Get_tag() == STOPTAG):
            print("worker {} stopping".format(comm.Get_rank()), flush=True)
            return
        # do work
        result = evolve_binary(input)
        # result = approximation_test(input)
        #print("output_file=",input.output_file)
        # output_file = open(input.output_file, 'w+')
        # print(result, file=output_file, flush=True)
        # indicate done with work by sending to something Master
        comm.send(result, dest=0)

########## Run the program
my_program()
