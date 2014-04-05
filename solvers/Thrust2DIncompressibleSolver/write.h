#include <fstream>
#include <iterator>
#include <string>


template <typename T>
void write_vector(  cusp::array1d<T, cusp::host_memory> V, 
                    const char* file_name){

    std::ofstream out(file_name);
    std::ostream_iterator <T> output_iterator(out, "\n");
    std::copy(V.begin(), V.end(), output_iterator);

}
