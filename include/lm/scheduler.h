/*
    Lightmetrica - Copyright (c) 2019 Hisanari Otsu
    Distributed under MIT license. See LICENSE file for details.
*/

#pragma once

#include "component.h"
#include <chrono>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)
LM_NAMESPACE_BEGIN(scheduler)

/*!
    \addtogroup scheduler
    @{
*/

/*!
    \brief Scheduler for rendering loop.

    \rst
    This interface provides an abstraction of the top-level rendering loop.
    \endrst
*/
class Scheduler : public Component {
public:
    /*!
        \brief Callback function for parallel loop.
        \param pixel_index Pixel index.
        \param sample_index Pixel sample index.
        \param threadid Thread index.
    */
    using ProcessFunc = std::function<void(long long pixel_index, long long sample_index, int threadid)>;

    /*!
        \brief Dispatch scheduler.
        \param process Callback function for parallel loop.
        \param before_func Callback function that will be called on each thread before any process callback is called
        \param after_func Callback function that will be called on thread after all process callbacks are called
        \return Processed samples per pixel.
    */
    virtual long long run(const ProcessFunc& process, const ProcessFunc& before_func = [&](auto,auto,auto) {} , const ProcessFunc& after_func = [&](auto,auto,auto) {}) const = 0;
};

/*!
    @}
*/

LM_NAMESPACE_END(scheduler)
LM_NAMESPACE_END(LM_NAMESPACE)
