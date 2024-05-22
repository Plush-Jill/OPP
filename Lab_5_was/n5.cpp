#include <thread>
#include <iostream>

void countNums1(long long begin, long long end, long long* counter) {
    long long num;

    for (long long i = begin; i <= end; ++i){
        if (/*begin == 1000000000 && */i % 100000000 == 0) {
            std::cout << std::this_thread::get_id() << ": " << i << std::endl;
        }
        num = i;
        long long sum1 = 0;
        long long sum2 = 0;
        long long tmp = 0;
        for (int k = 0; k < 10; ++k){
            tmp = num % 10;
            if (tmp % 2 == 0){
                sum2 += tmp;
            } else {
                sum1 += tmp;
            }
            num /= 10;
        }
        std::string string;
        if (sum1 < sum2) {
            string = std::to_string(sum1) + std::to_string(sum2);
        } else {
            string = std::to_string(sum2) + std::to_string(sum1);
        }
        //std::cout << string << std::endl;
        if (string == "3843") {
            ++(*counter);
        }

    }
}
void countNums2(long long begin, long long end, long long* counter) {
    long long num;

    for (long long i = begin; i <= end; ++i){
        if (/*begin == 1000000000 && */i % 100000000 == 0) {
            std::cout << std::this_thread::get_id() << ": " << i << std::endl;
        }
        num = i;
        long long sum1 = 0;
        long long sum2 = 0;
        long long tmp = 0;
        for (int k = 0; k < 10; ++k){
            tmp = num % 10;
            if (tmp % 2 == 0){
                sum2 += tmp;
            } else {
                sum1 += tmp;
            }
            num /= 10;
        }
        std::string string;
        if (sum1 < sum2) {
            string = std::to_string(sum1) + std::to_string(sum2);
        } else {
            string = std::to_string(sum2) + std::to_string(sum1);
        }
        if (string == "3843") {
            ++(*counter);
        }

    }
}
void countNums3(long long begin, long long end, long long* counter) {
    long long num;

    for (long long i = begin; i <= end; ++i){
        if (/*begin == 1000000000 && */i % 100000000 == 0) {
            std::cout << std::this_thread::get_id() << ": " << i << std::endl;
        }
        num = i;
        long long sum1 = 0;
        long long sum2 = 0;
        long long tmp = 0;
        for (int k = 0; k < 10; ++k){
            tmp = num % 10;
            if (tmp % 2 == 0){
                sum2 += tmp;
            } else {
                sum1 += tmp;
            }
            num /= 10;
        }
        std::string string;
        if (sum1 < sum2) {
            string = std::to_string(sum1) + std::to_string(sum2);
        } else {
            string = std::to_string(sum2) + std::to_string(sum1);
        }
        if (string == "3843") {
            ++(*counter);
        }

    }
}
void countNums4(long long begin, long long end, long long* counter) {
    long long num;

    for (long long i = begin; i <= end; ++i){
        if (/*begin == 1000000000 && */i % 100000000 == 0) {
            std::cout << std::this_thread::get_id() << ": " << i << std::endl;
        }
        num = i;
        long long sum1 = 0;
        long long sum2 = 0;
        long long tmp = 0;
        for (int k = 0; k < 10; ++k){
            tmp = num % 10;
            if (tmp % 2 == 0){
                sum2 += tmp;
            } else {
                sum1 += tmp;
            }
            num /= 10;
        }
        std::string string;
        if (sum1 < sum2) {
            string = std::to_string(sum1) + std::to_string(sum2);
        } else {
            string = std::to_string(sum2) + std::to_string(sum1);
        }
        if (string == "3843") {
            ++(*counter);
        }

    }
}
void countNums5(long long begin, long long end, long long* counter) {
    long long num;

    for (long long i = begin; i <= end; ++i){
        if (/*begin == 1000000000 && */i % 100000000 == 0) {
            std::cout << std::this_thread::get_id() << ": " << i << std::endl;
        }
        num = i;
        long long sum1 = 0;
        long long sum2 = 0;
        long long tmp = 0;
        for (int k = 0; k < 10; ++k){
            tmp = num % 10;
            if (tmp % 2 == 0){
                sum2 += tmp;
            } else {
                sum1 += tmp;
            }
            num /= 10;
        }
        std::string string;
        if (sum1 < sum2) {
            string = std::to_string(sum1) + std::to_string(sum2);
        } else {
            string = std::to_string(sum2) + std::to_string(sum1);
        }
        if (string == "3843") {
            ++(*counter);
        }

    }
}
void countNums6(long long begin, long long end, long long* counter) {
    long long num;

    for (long long i = begin; i <= end; ++i){
        if (/*begin == 1000000000 && */i % 100000000 == 0) {
            std::cout << std::this_thread::get_id() << ": " << i << std::endl;
        }
        num = i;
        long long sum1 = 0;
        long long sum2 = 0;
        long long tmp = 0;
        for (int k = 0; k < 10; ++k){
            tmp = num % 10;
            if (tmp % 2 == 0){
                sum2 += tmp;
            } else {
                sum1 += tmp;
            }
            num /= 10;
        }
        std::string string;
        if (sum1 < sum2) {
            string = std::to_string(sum1) + std::to_string(sum2);
        } else {
            string = std::to_string(sum2) + std::to_string(sum1);
        }
        if (string == "3843") {
            ++(*counter);
        }

    }
}

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    long long counter1 = 0;
    long long counter2 = 0;
    long long counter3 = 0;
    long long counter4 = 0;
    long long counter5 = 0;
    long long counter6 = 0;

    std::thread thread1(countNums1, 1000000000, 1700000000, &counter1);
    std::thread thread2(countNums1, 1700000000 + 1, 3400000000, &counter2);
    std::thread thread3(countNums1, 3400000000 + 1, 5000000000, &counter3);
    std::thread thread4(countNums1, 5000000000 + 1, 6700000000, &counter4);
    std::thread thread5(countNums1, 6700000000 + 1, 8400000000, &counter5);
    std::thread thread6(countNums1, 8400000000 + 1, 9999999999, &counter6);

    thread1.join();
    thread2.join();
    thread3.join();
    thread4.join();
    thread5.join();
    thread6.join();

    long long sum = counter1 + counter2 + counter3 + counter4 + counter5 + counter6;
    std::cout << sum << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}