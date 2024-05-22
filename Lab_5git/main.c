#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
#include <pthread.h>
#include <stdbool.h>


#define SUCCESS 0
#define ERROR   (-1)
#define TASK_COUNT              10000
#define TOTAL_SUM_WEIGHT        500000
#define REQUEST_TAG             0
#define RESPONSE_TAG            1
#define EMPTY_QUEUE_RESPONSE    (-1)
#define TERMINATION_SIGNAL      (-2)

typedef struct {
    int ID;
    int processID;
    int weight;
} Task;

typedef struct {
    Task *taskList;
    int capacity;
    int size;
    int popIndex;
} TaskQueue;

TaskQueue *createTaskQueue(int capacity) {
    TaskQueue *task_queue = malloc(sizeof(TaskQueue));

    if (task_queue == NULL) {
        return NULL;
    }

    Task *task_list = malloc (sizeof(Task) * capacity);
    if (task_list == NULL) {
        return NULL;
    }

    task_queue->taskList = task_list;
    task_queue->capacity = capacity;
    task_queue->size = 0;
    task_queue->popIndex = 0;

    return task_queue;
}
bool isTaskQueueEmpty(const TaskQueue *task_queue) {
    return task_queue->size == 0;
}
bool isTaskQueueFull(const TaskQueue *task_queue) {
    return task_queue->size == task_queue->capacity;
}
int pushIntoTaskQueue(TaskQueue *task_queue, Task task) {
    if (task_queue == NULL) {
        return ERROR;
    }

    if (isTaskQueueFull(task_queue)) {
        return ERROR;
    }

    int pushIndex = (task_queue->popIndex + task_queue->size) % task_queue->capacity;
    task_queue->taskList[pushIndex] = task;
    task_queue->size++;

    return SUCCESS;
}
int popFromTaskQueue(TaskQueue *task_queue, Task *task) {
    if (task_queue == NULL) {
        return ERROR;
    }

    if (isTaskQueueEmpty(task_queue)) {
        return ERROR;
    }

    *task = task_queue->taskList[task_queue->popIndex];
    task_queue->popIndex = (task_queue->popIndex + 1) % task_queue->capacity;
    task_queue->size--;

    return SUCCESS;
}
void destroyTaskQueue(TaskQueue **task_queue) {
    if (*task_queue == NULL) {
        return;
    }

    if ((*task_queue)->taskList == NULL) {
        return;
    }

    free((*task_queue)->taskList);
    free(*task_queue);

    *task_queue = NULL;
}


static int processID;
static int processCount;
static int startWeightSum = 0;
static int endWeightSum = 0;
bool termination = false;
TaskQueue* taskQueue_;

pthread_mutex_t mutex;
pthread_cond_t workerCondition;
pthread_cond_t receiverCondition;
static double globalRes = 0;

static inline void initTasks() {
    int minWeight = 2 * TOTAL_SUM_WEIGHT / (TASK_COUNT * (processCount + 1));
    int taskID = 1;

    for (int i = 0; i < TASK_COUNT; ++i) {
        Task task = {
                .ID = taskID,
                .processID = processID,
                .weight = minWeight * (i % processCount + 1)
        };

        if (i % processCount == processID) {
            pushIntoTaskQueue(taskQueue_, task);
            ++taskID;
            startWeightSum += task.weight;
        }
    }
}
static inline void executeCurrentTasks() {
    while (true) {
        Task task;

        pthread_mutex_lock(&mutex);
        if (isTaskQueueEmpty(taskQueue_)) {
            pthread_mutex_unlock(&mutex);
            break;
        }
        popFromTaskQueue(taskQueue_, &task);
        pthread_mutex_unlock(&mutex);

        for (int i = 0; i < task.weight; ++i) {
            for (int j = 0; j < 20000; ++j) {
                globalRes += sqrt(sqrt(sqrt(sqrt(i))));
            }
        }

        endWeightSum += task.weight;
    }
}
void *workerStart() {
    initTasks();

    MPI_Barrier(MPI_COMM_WORLD);

    while (true) {
        executeCurrentTasks();
        pthread_mutex_lock(&mutex);
        while (isTaskQueueEmpty(taskQueue_) && !termination) {
            pthread_cond_signal(&receiverCondition);
            pthread_cond_wait(&workerCondition, &mutex);
        }

        if (termination) {
            pthread_mutex_unlock(&mutex);
            break;
        }
        pthread_mutex_unlock(&mutex);
    }

    printf("Worker %d finished\n", processID);
    pthread_exit(NULL);
}
void *receiverStart() {
    int terminationSignal = TERMINATION_SIGNAL;

    while (!termination) {
        int receivedTasks = 0;
        Task task;

        pthread_mutex_lock(&mutex);
        while (!isTaskQueueEmpty(taskQueue_)) {
            pthread_cond_wait(&receiverCondition, &mutex);
        }
        pthread_mutex_unlock(&mutex);

        for (int i = 0; i < processCount; ++i) {
            if (i == processID) {
                continue;
            }

            MPI_Send(&processID, 1, MPI_INT, i, REQUEST_TAG, MPI_COMM_WORLD);
            MPI_Recv(&task, sizeof(task), MPI_BYTE, i , RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (task.ID != EMPTY_QUEUE_RESPONSE) {
                pthread_mutex_lock(&mutex);
                pushIntoTaskQueue(taskQueue_, task);
                pthread_mutex_unlock(&mutex);

                ++receivedTasks;
            }
        }

        if (receivedTasks == 0) {
            pthread_mutex_lock(&mutex);
            termination = true;
            pthread_mutex_unlock(&mutex);
        }

        pthread_mutex_lock(&mutex);
        pthread_cond_signal(&workerCondition);
        pthread_mutex_unlock(&mutex);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Send(&terminationSignal, 1, MPI_INT, processID, REQUEST_TAG, MPI_COMM_WORLD);
    pthread_exit(NULL);
}
void *senderStart() {
    while (true) {
        int receiverProcessID;

        Task task;

        MPI_Recv(&receiverProcessID, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (receiverProcessID == TERMINATION_SIGNAL) {
            break;
        }

        pthread_mutex_lock(&mutex);
        if (!isTaskQueueEmpty(taskQueue_)) {
            popFromTaskQueue(taskQueue_, &task);
        } else {
            task.ID = EMPTY_QUEUE_RESPONSE;
            task.weight = 0;
            task.processID = processID;
        }
        pthread_mutex_unlock(&mutex);

        MPI_Send(&task, sizeof(task), MPI_BYTE, receiverProcessID, RESPONSE_TAG, MPI_COMM_WORLD);
    }

    pthread_exit(NULL);
}

int main(int argc, char **argv) {
    int required = MPI_THREAD_MULTIPLE;
    int provided;
    double beginningTime;
    double endingTime;
    pthread_t workerThread;
    pthread_t receiverThread;
    pthread_t senderThread;

    MPI_Init_thread(&argc, &argv, required, &provided);
    if (provided != required) {
        return EXIT_FAILURE;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    taskQueue_ = createTaskQueue(TASK_COUNT);

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&workerCondition, NULL);
    pthread_cond_init(&receiverCondition, NULL);

    beginningTime = MPI_Wtime();
    pthread_create(&workerThread, NULL, workerStart, NULL);
    pthread_create(&receiverThread, NULL, receiverStart, NULL);
    pthread_create(&senderThread, NULL, senderStart, NULL);

    pthread_join(workerThread, NULL);
    pthread_join(receiverThread, NULL);
    pthread_join(senderThread, NULL);
    endingTime = MPI_Wtime();

    double time = endingTime - beginningTime;
    double finalTime = 0;
    MPI_Reduce(&time, &finalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Summary weight %d - start: %d, actual: %d\n", processID, startWeightSum, endWeightSum);
    MPI_Barrier(MPI_COMM_WORLD);
    if (processID == 0) {
        printf("Time: %lf\n", finalTime);
    }

    destroyTaskQueue(&taskQueue_);
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&workerCondition);
    pthread_cond_destroy(&receiverCondition);
    MPI_Finalize();

    return EXIT_SUCCESS;
}