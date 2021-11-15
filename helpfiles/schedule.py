from mpi4py import MPI
import numpy as np
import time

def genTasks(numTasks):
    np.random.seed(1000)  # run the same set of timed tasks
    return np.random.randint(low=1, high=20, size=numTasks)

# tags that can be applied to messages
WORKTAG = 1
STOPTAG = 2

def my_program():
    comm = MPI.COMM_WORLD
    id = comm.Get_rank()            #number of the process
    numProc = comm.Get_size()  #total number of processes

    if (id == 0) :
        numTasks = (numProc-1)*4 # avg 4 tasks per worker process
        workTimes = genTasks(numTasks)
        print(workTimes, flush=True)
        superwise_work(workTimes, comm, numProc)
    else:
        do_work(comm)

def superwise_work(workTimes, comm, numProc):
    totalWork = workTimes.size
    workcount = 0
    recvcount = 0
    # send out the first tasks to all workers
    for id in range(1, numProc):
        if workcount < totalWork:
            work=workTimes[workcount]
            comm.send(work, dest=id, tag=WORKTAG)
            workcount += 1
            print("master sent {} to {}".format(work, id), flush=True)

    # while there is still work,
    # receive result/requests for more work from a worker
    while (workcount < totalWork) :
        # receive next finished result
        stat = MPI.Status()
        workTime = comm.recv(source=MPI.ANY_SOURCE, status=stat)
        recvcount += 1
        workerId = stat.Get_source()
        print("master received {} from {}".format(workTime, workerId), flush=True)
        #send next work
        comm.send(workTimes[workcount], dest=workerId, tag=WORKTAG)
        workcount += 1
        print("master sent {} to {}".format(work, workerId), flush=True)

    # Receive results for outstanding work requests.
    while (recvcount < totalWork):
        stat = MPI.Status()
        workTime = comm.recv(source=MPI.ANY_SOURCE, status=stat)
        recvcount += 1
        workerId = stat.Get_source()
        print("end: master received {} from {}".format(workTime, workerId), flush=True)

    # Tell all workers to stop
    for id in range(1, numProc):
        comm.send(-1, dest=id, tag=STOPTAG)


def do_work(comm):
    # loop till the tag says to stop
    while(True):
        stat = MPI.Status()
        waitTime = comm.recv(source=0, tag=MPI.ANY_TAG, status=stat)
        print("worker {} got {}".format(comm.Get_rank(), waitTime), flush=True)
        if (stat.Get_tag() == STOPTAG):
            print("worker {} stopping".format(comm.Get_rank()), flush=True)
            return
        # do work
        time.sleep(waitTime)
        # indicate done with work by sending to something Master
        comm.send(waitTime, dest=0)

########## Run the program
my_program()
