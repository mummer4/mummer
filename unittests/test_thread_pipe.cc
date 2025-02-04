#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <gtest/gtest.h>

#include <thread_pipe.hpp>

namespace {
    // Sizes written and block size
    static const size_t times = 100000;
    static const size_t size = 67;
    static const ssize_t block = 1024;

    // Thread function: Output
    void producer(thread_pipe::ostream_buffered* output, int thread_id) {
        auto it = output->begin();

        for(size_t i = 0; i < times; ++i) {
            *it << thread_id;
            for(size_t j = 0; j < size; ++j)
              *it << ' ' << i;
            *it << '\t';
            if(it->tellp() > block)
              ++it;
        }
        it.done();
    }

    TEST(ThreadPipe, MultipleProducers) {
        static const char* file = "multipleproducers";
        static const int nb_threads = 4;

        { // Write content to file
            std::ofstream os(file);
            thread_pipe::ostream_buffered output(os);

            std::vector<std::thread> threads;
            for(int i = 0; i < nb_threads; ++i)
                threads.push_back(std::thread(producer, &output, i));

            for(auto& th : threads)
                th.join();

            EXPECT_TRUE(os.good());
        }

        { // Read and check content in file. It is tab separated
            std::ifstream is(file);
            std::vector<std::string> content;
            std::string block;

            while(std::getline(is, block, '\t')) {
                content.push_back(block);
            }

            EXPECT_TRUE(is.eof());
            EXPECT_EQ(times * (size_t)nb_threads, content.size());

            // Expect every thread to have written lines in order, each line
            // containing "times" order value.
            std::vector<size_t> indices(nb_threads, 0);
            std::istringstream iss;
            int thid;
            size_t count;
            for(const auto& l : content) {
                iss.str(l);
                iss.clear();
                iss >> thid;
                EXPECT_TRUE(iss.good());
                EXPECT_GE(thid, 0);
                EXPECT_LT(thid, nb_threads);
                for(size_t i = 0; i < size; ++i) {
                    EXPECT_TRUE(iss.good()) << "thid " << this << " i " << i;
                    iss >> count;
                    EXPECT_EQ(indices[thid], count) << "thid " << thid;
                }
                iss >> count;
                EXPECT_TRUE(iss.eof());
                ++indices[thid];
            }

            for(const auto v : indices)
                EXPECT_EQ(times, v);
        }
    }

} // namespace
