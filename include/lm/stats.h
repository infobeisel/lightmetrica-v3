#pragma once
#include "component.h"
#include <unordered_map>
#include <functional>
#include <memory>
#include <lm/json.h>
#include <mutex>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)

namespace stats {

/*!
    \addtogroup stats
    \brief namespace containing functions to collect statistics (e.g. during a render) thread-locally and merge them to global statistics afterwards.
    
    typically one first creates a type that serves as Tag to identify a group of statistics, e.g.
    class TrianglesHit {};
    then one calls the clear functions to clear statistics data that has been accumulated so far (statistics are collected in a global manner!)
    clearGlobal<TrianglesHit,int,int>(); 
    and inside each thread
    clear<TrianglesHit,int,int>();
    In the example calls above, the statistics have the tag "TrianglesHit" and each data entry has an int as a key and an int as a value (could be triangle id and number of hits for example)
    
    Then one should use the reserve<...> methods to prepare the module for the amount of incoming data
    
    To actually collect statistics, one then first uses the thread local functions
    set<...>(), add<...>(), and the general insert<...>()
    and
    get<...>() to receive thread-locally collected data. 
    After the collection is finished, one needs to merge the thread local statistics to the global statistics calling 
    mergeToGlobal<...>().
    Then one can call getGlobal<...>() to receive the total statistics data and
    JSONGatherer::to_json(...) to dump all gathered statistics of any type to one json

    Alternatively, if the calls above are too slow (they all use hashing), one can accumulate
    data by calling the enqueue<...>() method. After the performance-critical work is done, one can flushQueue<...>(), giving a lambda that processes
    the accumulated data . For example, the lambda could contain some set<...>() or add<..>() calls, and after the flushQueue there could be a mergeToGlobal<...>() call again 

    Guide to get meaningful timings although collecting statistics when rendering:
    make sure to move the clear*(), reserve*(), mergeToGlobal() (and the flushQueue()) calls outside the time critical blocks and to only include set()/ add()/ insert()/ enqueue() calls.
    the parallel::foreach provides a before and after process function where you can set up timers for the above calls, which you can in the end subtract from the total time
    
    
    @{
*/

/*!
    \brief class to be able to gather jsons of all potential Statistics<...> instantiations
*/
class JSONGatherer {
    private:
        JSONGatherer() {}
        ~JSONGatherer(){}
        std::vector<std::function<void(lm::Json &)>> toJsonFs;
    public:
        void registerToJsonF(std::function<void(lm::Json &)> f) {
            toJsonFs.push_back(f);
        }
        void to_json(lm::Json & j) {
            lm::Json merged;
            lm::Json toMerge;
            
            for(auto f : toJsonFs) {
                f(toMerge);
                merged.push_back(toMerge);
            }

            j = merged;
        }
        static JSONGatherer & instance() {
            static JSONGatherer s;
            return s;
        }
};

/*!
    \brief stats internal implementation namespace
*/
namespace internals {
    


    template <typename Tag, typename Key, typename Value>
    /*!
    \brief internal statistics implementation,(one static and one thread_local) singleton, to access a std::unordered_map
    */
    class Statistics {
        private:
            Statistics() {}
            Statistics(bool registerJsonGather) {
                JSONGatherer::instance().registerToJsonF( [&](lm::Json & j) {
                    Statistics<Tag,Key,Value>::accessMainInstance([&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
                        j = mainInstance.data;
                    });
                });
            }
            ~Statistics() {}
            
        public:
            std::unordered_map<Key,Value> data;
            std::vector<Key> keyQueue;
            std::vector<Value> valueQueue;
            static void accessMainInstance(std::function<void(Statistics<Tag,Key,Value>&)> accessor) {
                static std::mutex access;
                access.lock();
                static Statistics<Tag, Key,Value> s(true);
                accessor(s);
                access.unlock();
            } 
            static Statistics<Tag,Key,Value> & threadInstance() {
                thread_local  Statistics<Tag,Key,Value> s;
                return s;
            }
            void shrinkClear() {
                data = std::unordered_map<Key,Value>();
                keyQueue = std::vector<Key>();
                valueQueue = std::vector<Value>();
            }

    };

    template <typename Tag, typename Value>
    /*!
    \brief internal statistics implementation,(one static and one thread_local) singleton, to access a std::unordered_map
    */
    class SlotStats {
        private:
            SlotStats() {}
            SlotStats(bool registerJsonGather) {
                JSONGatherer::instance().registerToJsonF( [&](lm::Json & j) {
                    SlotStats<Tag,Value>::accessMainInstance([&] (internals::SlotStats<Tag,Value> & mainInstance) {
                        j = mainInstance.data;
                    });
                });
            }
            ~SlotStats() {}


            
        public:
            std::vector<Value> data;
            static void accessMainInstance(std::function<void(SlotStats<Tag,Value>&)> accessor) {
                static std::mutex access;
                access.lock();
                accessor(mainInstance());
                access.unlock();
            } 

            static SlotStats<Tag,Value> & mainInstance() {
                static SlotStats<Tag,Value> s(true);
                return s;
            }
            
            static SlotStats<Tag,Value> & threadInstance() {
                thread_local  SlotStats<Tag,Value> s;
                return s;
            }
            void shrinkClear() {
                data = std::vector<Value>();
            }

    };
}



template <typename Tag, typename Key, typename Value>
/*!
    \brief parses statistics from json, using the given converter lambda to parse a json entry into a Key-Value pair
*/
LM_PUBLIC_API void from_json(lm::Json & j, std::function<std::pair<Key,Value>(lm::Json & keyJ, lm::Json & valJ)> converter) {
    internals::Statistics<Tag, Key,Value>::accessMainInstance( 
        [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
            for(auto & kv : j.items()) {
                for(nlohmann::basic_json<>::iterator it = kv.value().begin(); it != kv.value().end(); ++it) {
                    auto p = converter(it.value()[0],it.value()[1]);
                    mainInstance.data[kv.key()][p.first] = p.second;
                }
            }
    });
}




template <typename Tag, typename Key, typename Value>
/*!
    \brief merges this thread's locally accumulated statistics into the global statistics, using the given lambda (which receives the current global Value and the thread local Value to merge in, should return the merged Value)
*/
LM_PUBLIC_API void mergeToGlobal(std::function<Value(Value &, Value &)> merger ) {
    internals::Statistics<Tag, Key,Value>::accessMainInstance( 
        [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
        for( auto entry : internals::Statistics<Tag, Key,Value>::threadInstance().data) {
            mainInstance.data[entry.first] = merger(mainInstance.data[entry.first] , entry.second);
        }
    });
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief merges this thread's locally accumulated statistics into the global statistics, using the given lambda (which receives the current global Value and the thread local Value to merge in, should return reference to the merged Value)
*/
LM_PUBLIC_API void mergeToGlobal(std::function<Value&(Value &, Value &)> merger ) {
    internals::Statistics<Tag, Key,Value>::accessMainInstance( 
        [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
        for( auto entry : internals::Statistics<Tag, Key,Value>::threadInstance().data) {
            mainInstance.data[entry.first] = merger(mainInstance.data[entry.first] , entry.second);
        }
    });
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief receives global statistics that have been accumulated so far
*/
LM_PUBLIC_API std::unordered_map<Key,Value> getGlobal() {
    std::unordered_map<Key,Value> ret;
    internals::Statistics<Tag, Key,Value>::accessMainInstance(
         [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
        ret = mainInstance.data;
    });
    return ret;
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief receives global statistics that have been accumulated so far
*/
LM_PUBLIC_API std::unordered_map<Key,Value> & getGlobalRef() {
    std::unordered_map<Key,Value> * ret;
    internals::Statistics<Tag, Key,Value>::accessMainInstance(
         [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
        ret = &mainInstance.data;
    });
    return *ret;
}






template <typename Tag, typename Key, typename Value>
/*!
    \brief receives global statistics that have been accumulated so far, belonging to given key
*/
LM_PUBLIC_API Value getGlobal( Key key) {
    Value ret;
    internals::Statistics<Tag, Key,Value>::accessMainInstance(
         [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
        ret = mainInstance.data[key];
    });
    return ret;
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief receives global statistics that have been accumulated so far, accessed in a not threadsafe manner
*/
LM_PUBLIC_API std::unordered_map<Key,Value> & getGlobalRefUnsafe(Key key) {
    return internals::Statistics<Tag, Key,Value>::mainInstance().data[key];
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief clears all global statistics
*/
LM_PUBLIC_API void clearGlobal() {
    internals::Statistics<Tag, Key,Value>::accessMainInstance(
         [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
            mainInstance.shrinkClear();
    });
}

template <typename Tag, typename Key, typename Value>
/*!
    \brief preallocates space for given entrycount
*/
LM_PUBLIC_API void reserveGlobal(size_t entryCount) {
    internals::Statistics<Tag, Key,Value>::accessMainInstance(
         [&] (internals::Statistics<Tag, Key,Value> & mainInstance) {
         mainInstance.data.reserve(entryCount);
    });
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief preallocates space for given entrycount, thread local
*/
LM_PUBLIC_API void reserve(size_t entryCount) {
    internals::Statistics<Tag, Key,Value>::threadInstance().data.reserve(entryCount);
    internals::Statistics<Tag, Key,Value>::threadInstance().keyQueue.reserve(entryCount);
    internals::Statistics<Tag, Key,Value>::threadInstance().valueQueue.reserve(entryCount);
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief sets a thread local statistics data point, overwrites previous value with the same key
*/
LM_PUBLIC_API void set( Key key, Value val) {
    internals::Statistics<Tag, Key,Value>::threadInstance().data[key] = val;
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief inserts a thread local statistics data point, given key, using the given inserter (which receives the current value and the new given value and should return value to insert)
*/
LM_PUBLIC_API void insert(Key key, Value value, std::function<Value(Value &, Value &)> inserter) {
    internals::Statistics<Tag, Key,Value>::threadInstance().data[key] = 
    inserter(
        internals::Statistics<Tag, Key,Value>::threadInstance().data[key],
        value);
}

template <typename Tag, typename Key, typename Value>
/*!
    \brief updates a thread local statistics data point, given key, using the given updater (which receives the current value (will be default constructed if not yet existent) and should return the reference to the updated original value)
*/
LM_PUBLIC_API void update(Key key, std::function<void(Value &)> updater) {
    updater(internals::Statistics<Tag, Key,Value>::threadInstance().data[key]);
}





template <typename Tag, typename Key, typename Value>
/*!
    \brief receives a thread local statistics data point, given key
*/
LM_PUBLIC_API Value get( Key & key) {
    return internals::Statistics<Tag, Key,Value>::threadInstance().data[key];
}

template <typename Tag, typename Key, typename Value>
/*!
    \brief receives a reference to a thread local statistics data point , given key
*/
LM_PUBLIC_API Value & getRef( Key key) {
    return internals::Statistics<Tag, Key,Value>::threadInstance().data[key];
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief receives a reference to thread local statistics 
*/
LM_PUBLIC_API std::unordered_map<Key,Value> & getRef() {
    return internals::Statistics<Tag, Key,Value>::threadInstance().data;
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief receives a copy to thread local statistics 
*/
LM_PUBLIC_API std::unordered_map<Key,Value> get() {
    return internals::Statistics<Tag, Key,Value>::threadInstance().data;
}



template <typename Tag, typename Key, typename Value,
        typename std::enable_if< std::is_arithmetic< Value >::value,Value >::type * = nullptr>
/*!
    \brief adds a given value to a previously given accumulated value, thread local
*/
LM_PUBLIC_API Value add( Key key, Value val) {
    internals::Statistics<Tag, Key,Value>::threadInstance().data[key] += val;
}
template <typename Tag, typename Key, typename Value>
/*!
    \brief erases a statistics entry (thread local)
*/
LM_PUBLIC_API void erase( Key key) {
    internals::Statistics<Tag, Key,Value>::threadInstance().data.erase(key);
}

template <typename Tag, typename Key, typename Value>
/*!
    \brief clears all statistics (thread local)
*/
LM_PUBLIC_API void clear() {
    internals::Statistics<Tag, Key,Value>::threadInstance().shrinkClear();
}

template <typename Tag, typename Key, typename Value>
/*!
    \brief enqueues the given data to the unprocessed statistics queue (thread local). these entries are being ignored by get* calls (use flushQueue() to actually include them into the statistics )
*/
LM_PUBLIC_API void enqueue(Key && k, Value && v) {
    internals::Statistics<Tag, Key,Value>::threadInstance().keyQueue.emplace_back(k);
    internals::Statistics<Tag, Key,Value>::threadInstance().valueQueue.emplace_back(v);
}


template <typename Tag, typename Key, typename Value>
/*!
    \brief flushes the unprocessed statistics queue, processing it with the given lambda(thread local)
*/
LM_PUBLIC_API void flushQueue( std::function<void(Key & k, Value & v)> processor) {

    auto keyIt = internals::Statistics<Tag, Key,Value>::threadInstance().keyQueue.begin();
    auto valIt = internals::Statistics<Tag, Key,Value>::threadInstance().valueQueue.begin();
    auto kEnd = internals::Statistics<Tag, Key,Value>::threadInstance().keyQueue.end();
    auto vEnd = internals::Statistics<Tag, Key,Value>::threadInstance().valueQueue.end();
    while(keyIt != kEnd && valIt != vEnd) {
        processor(*keyIt,*valIt);
        keyIt++;
        valIt++;    
    }
    internals::Statistics<Tag, Key,Value>::threadInstance().keyQueue.clear();
    internals::Statistics<Tag, Key,Value>::threadInstance().valueQueue.clear();
}



template <typename Tag, typename Value>
/*!
    \brief sets a thread local statistics data point, overwrites previous value with the same index
*/
LM_PUBLIC_API void set( size_t index, Value val) {
    if(internals::SlotStats<Tag,Value>::threadInstance().data.size() < index + 1)
        internals::SlotStats<Tag,Value>::threadInstance().data.resize(index + 1);
    internals::SlotStats<Tag,Value>::threadInstance().data[index] = val;
}
 
template <typename Tag, typename Value>
/*!
    \brief receives a thread local statistics data point, given index
*/
LM_PUBLIC_API Value get( size_t index) {
    if(internals::SlotStats<Tag,Value>::threadInstance().data.size() < index + 1)
        internals::SlotStats<Tag,Value>::threadInstance().data.resize(index + 1);
    return internals::SlotStats<Tag,Value>::threadInstance().data[index];
}

template <typename Tag, typename Value>
/*!
    \brief receives a reference to a thread local statistics data point , given index
*/
LM_PUBLIC_API Value & getRef(size_t index) {
    if(internals::SlotStats<Tag,Value>::threadInstance().data.size() < index + 1)
        internals::SlotStats<Tag,Value>::threadInstance().data.resize(index + 1);
    return internals::SlotStats<Tag,Value>::threadInstance().data[index];
}


template <typename Tag, typename Value>
/*!
    \brief preallocates space for given entrycount, thread local
*/
LM_PUBLIC_API void reserve(size_t entryCount) {
    internals::SlotStats<Tag, Value>::threadInstance().data.resize(entryCount);
}

/*!
    @}
*/
}
LM_NAMESPACE_END(LM_NAMESPACE)