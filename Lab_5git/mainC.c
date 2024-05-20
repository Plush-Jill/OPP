#include <mpi/mpi.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#define FNORM       "\x1B[0m"
#define FBLACK      "\x1B[30m"
#define FRED        "\x1B[31m"
#define FGREEN      "\x1B[32m"
#define FYELLOW     "\x1B[33m"
#define FBLUE       "\x1B[34m"
#define FMAGENTA    "\x1B[35m"
#define FCYAN       "\x1B[36m"
#define FWHITE      "\x1B[37m"

#define TASK_COUNT              2000
#define REQUEST_TAG             0
#define RESPONSE_TAG            1
#define EMPTY_QUEUE_RESPONSE    (-1)
#define TERMINATION_SIGNAL      (-2)
#define SUCCESS 0
#define ERROR   (-1)

struct task_t {
    int ID;
    int processID;
    int weight;
};
struct taskQueue_t {
    struct task_t* data;
    int capacity;
    int count;
    int popIndex;
};

int processCount;
int processID;
int procSumWeight = 0;
bool termination = false;
struct taskQueue_t* taskQueue;

pthread_mutex_t mutex;
pthread_cond_t workerCondition;
pthread_cond_t receiverCondition;

void* workerStart(void* args);
void* receiverStart(void* args);
void* senderStart(void* args);

struct taskQueue_t* createTaskQueue(int capacity);
bool isTaskQueueEmpty(const struct taskQueue_t* queue);
bool isTaskQueueFull(const struct taskQueue_t* queue);
int taskQueuePush(struct taskQueue_t* queue, struct task_t task);
int taskQueuePop(struct taskQueue_t* queue, struct task_t* task);
void taskQueueDestroy(struct taskQueue_t* *queue);


int main(int argc, char* argv[]) {
    int required = MPI_THREAD_MULTIPLE;
    int provided;
    double beginningTime;
    double endingTime;
    pthread_t workerThread;
    pthread_t receiverThread;
    pthread_t senderThread;

    // Initialize MPI environment
    MPI_Init_thread(&argc, &argv, required, &provided);
    if (provided != required) {
        return EXIT_FAILURE;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &processID);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    // Create task queue
    taskQueue = createTaskQueue(TASK_COUNT);
    fprintf(stdout, "%d %d %d\n",
            taskQueue->capacity, taskQueue->count, taskQueue->popIndex);
    // Initialize mutex and condition variable
    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&workerCondition, NULL);
    pthread_cond_init(&receiverCondition, NULL);

    // Start worker and sender thread
    beginningTime = MPI_Wtime();
    pthread_create(&workerThread, NULL, workerStart, NULL);
    pthread_create(&receiverThread, NULL, receiverStart, NULL);
    pthread_create(&senderThread, NULL, senderStart, NULL);

    pthread_join(workerThread, NULL);
    pthread_join(receiverThread, NULL);
    pthread_join(senderThread, NULL);
    endingTime = MPI_Wtime();

    // Print result
    MPI_Barrier(MPI_COMM_WORLD);
    printf(FGREEN"Summary weight %d: %lf\n" FNORM, processID, procSumWeight * 1E-6);
    MPI_Barrier(MPI_COMM_WORLD);
    if (processID == 0) {
        printf(FGREEN"Time: %lf\n" FNORM, endingTime - beginningTime);
    }

    // Clear resources
    taskQueueDestroy(&taskQueue);
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&workerCondition);
    pthread_cond_destroy(&receiverCondition);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

static inline void initTasks() {
    // Total sum of task weights does not change
    // For each process, tasks have a weight: min_weight * (process_id + 1)
    // min_weight = (TOTAL_SUM_WEIGHT * n) / (TASK_COUNT * S_n)
    // n - process count
    // S_n - sum of an arithmetic progression 1, 2, ... , n

    const int TOTAL_SUM_WEIGHT = 50000000;
    int minWeight = 2 * TOTAL_SUM_WEIGHT / (TASK_COUNT * (processCount + 1));
    int taskID = 1;

    for (int i = 0; i < TASK_COUNT; ++i) {
        // Create task
        struct task_t task = {
                .ID = taskID,
                .processID = processID,
                .weight = minWeight * (i % processCount + 1)
        };

        if (i % processCount == processID) {
            taskQueuePush(taskQueue, task);
            ++taskID;
            procSumWeight += task.weight;
        }
    }
}
static inline void executeTasks() {
    while (true) {
        struct task_t task;

        pthread_mutex_lock(&mutex);
        if (isTaskQueueEmpty(taskQueue)) {
            pthread_mutex_unlock(&mutex);
            break;
        }
        taskQueuePop(taskQueue, &task);
        pthread_mutex_unlock(&mutex);

        printf(FBLUE"Worker %d executing task %d of process %d and weight %d\n"FNORM,
               processID,
               task.ID,
               task.processID,
               task.weight);
        usleep(task.weight);
    }
}
void* workerStart(void* args) {
    initTasks();

    // Worker start synchronization
    MPI_Barrier(MPI_COMM_WORLD);

    while (true) {
        executeTasks();

        pthread_mutex_lock(&mutex);
        if (isTaskQueueEmpty(taskQueue) && !termination){
            pthread_cond_signal(&receiverCondition);
            pthread_cond_wait(&workerCondition, &mutex);
        }

        if (termination) {
            pthread_mutex_unlock(&mutex);
            break;
        }
        pthread_mutex_unlock(&mutex);
    }

    printf(FBLUE"Worker %d finished\n"FNORM, processID);
    pthread_exit(NULL);
}
void* receiverStart(void* args) {
    int terminationSignal = TERMINATION_SIGNAL;

    while (!termination) {
        int receivedTasks = 0;
        struct task_t task;

        pthread_mutex_lock(&mutex);
        if (!isTaskQueueEmpty(taskQueue)) {
            pthread_cond_wait(&receiverCondition, &mutex);
        }
        pthread_mutex_unlock(&mutex);

        for (int i = 0; i < processCount; ++i) {
            if (i == processID) {
                continue;
            }

            printf(FYELLOW"Receiver %d sent request to process %d\n"FNORM, processID, i);
            MPI_Send(&processID, 1, MPI_INT, i, REQUEST_TAG, MPI_COMM_WORLD);
            MPI_Recv(&task, sizeof(task), MPI_BYTE, i, RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (task.ID != EMPTY_QUEUE_RESPONSE) {
                printf(FYELLOW"Receiver %d received task %d from process %d\n"FNORM, processID, task.ID, i);

                pthread_mutex_lock(&mutex);
                taskQueuePush(taskQueue, task);
                pthread_mutex_unlock(&mutex);

                ++receivedTasks;
            } else {
                printf(FYELLOW"Receiver %d received empty queue response from process %d\n"FNORM, processID, i);
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

    // Receiver destruction synchronization
    MPI_Barrier(MPI_COMM_WORLD);

    printf(FYELLOW"Receiver %d sent termination signal\n"FNORM, processID);
    MPI_Send(&terminationSignal, 1, MPI_INT, processID, REQUEST_TAG, MPI_COMM_WORLD);

    printf(FYELLOW"Receiver %d finished\n"FNORM, processID);
    pthread_exit(NULL);
}
void* senderStart(void* args) {
    while (true) {
        int receiveProcessID;
        struct task_t task;

        printf(FMAGENTA"Sender %d waiting for request\n"FNORM, processID);
        MPI_Recv(&receiveProcessID,1,MPI_INT,MPI_ANY_SOURCE,REQUEST_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        if (receiveProcessID == TERMINATION_SIGNAL){
            printf(FMAGENTA"Sender %d received termination signal\n"FNORM, processID);
            break;
        }

        printf(FMAGENTA"Sender %d received request from process %d\n"FNORM, processID, receiveProcessID);

        pthread_mutex_lock(&mutex);
        if (!isTaskQueueEmpty(taskQueue)) {
            taskQueuePop(taskQueue, &task);
            printf(FMAGENTA"Sender %d sent task %d of process %d to process %d\n"FNORM, processID, task.ID, task.processID, receiveProcessID);
        } else {
            task.ID = EMPTY_QUEUE_RESPONSE;
            task.weight = 0;
            task.processID = processID;
            printf(FMAGENTA"Sender %d sent empty queue response to process %d\n"FNORM, processID, receiveProcessID);
        }
        pthread_mutex_unlock(&mutex);

        MPI_Send(&task, sizeof(task), MPI_BYTE, receiveProcessID, RESPONSE_TAG, MPI_COMM_WORLD);
    }

    printf(FMAGENTA"Sender %d finished\n"FNORM, processID);
    pthread_exit(NULL);
}


struct taskQueue_t* createTaskQueue(int capacity) {
    struct taskQueue_t* queue = malloc(sizeof(struct taskQueue_t));
    if (queue == NULL) {
        return NULL;
    }

    struct task_t* data = malloc(sizeof(struct task_t) * capacity);
    if (data == NULL) {
        return NULL;
    }

    queue->data = data;
    queue->capacity = capacity;
    queue->count = 0;
    queue->popIndex = 0;

    return queue;
}
bool isTaskQueueEmpty(const struct taskQueue_t* queue) {
    return queue->count == 0;
}
bool isTaskQueueFull(const struct taskQueue_t* queue) {
    return queue->count == queue->capacity;
}
int taskQueuePush(struct taskQueue_t* queue, struct task_t task) {
    if (queue == NULL) {
        return ERROR;
    }

    if (isTaskQueueFull(queue)) {
        return ERROR;
    }

    int push_index = (queue->popIndex + queue->count) % queue->capacity;
    queue->data[push_index] = task;
    ++queue->count;

    return SUCCESS;
}
int taskQueuePop(struct taskQueue_t* queue, struct task_t* task) {
    if (queue == NULL) {
        return ERROR;
    }

    if (isTaskQueueEmpty(queue)) {
        return ERROR;
    }

    *task = queue->data[queue->popIndex];
    queue->popIndex = (queue->popIndex + 1) % queue->capacity;
    queue->count--;

    return SUCCESS;
}
void taskQueueDestroy(struct taskQueue_t* *queue) {
    if (*queue == NULL) {
        return;
    }

    if ((*queue)->data == NULL) {
        return;
    }

    free((*queue)->data);
    free(*queue);

    * queue = NULL;
}