/**
******************************************************************************
* @file:   barrier.hpp
*
* @author: Marvin Williams
* @date:   2021/10/28 14:26
* @brief:
*******************************************************************************
**/

#pragma once
#ifndef BARRIER_HPP_INCLUDED
#define BARRIER_HPP_INCLUDED

#include <pthread.h>

#include <exception>
#include <iostream>
#include <system_error>

namespace utils {

class barrier {
    pthread_barrier_t b_{};

   public:
    explicit barrier(int num_threads) {
        assert(num_threads > 0);
        if (int rc = pthread_barrier_init(&b_, nullptr, static_cast<unsigned int>(num_threads)); rc != 0) {
            throw std::system_error{rc, std::system_category(), "Failed to create barrier: "};
        }
    }

    barrier(barrier const &) = delete;
    barrier &operator=(barrier const &) = delete;
    barrier(barrier &&) = delete;
    barrier &operator=(barrier &&) = delete;

    void wait() {
        if (int rc = pthread_barrier_wait(&b_); rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
            throw std::system_error{rc, std::system_category(), "Failed to wait for barrier: "};
        }
    }

    template <typename CompletionFunc>
    void wait(CompletionFunc f) {
        if (int rc = pthread_barrier_wait(&b_); rc == PTHREAD_BARRIER_SERIAL_THREAD) {
            f();
        } else if (rc != 0) {
            throw std::system_error{rc, std::system_category(), "Failed to wait for barrier: "};
        }
    }

    ~barrier() noexcept {
        if (int rc = pthread_barrier_destroy(&b_); rc != 0) {
            std::cerr << "Failed to destroy barrier: " << std::system_category().message(rc) << " (Error " << rc
                      << ")\n";
            std::terminate();
        }
    }
};

}  // namespace utils

#endif  //! BARRIER_HPP_INCLUDED
