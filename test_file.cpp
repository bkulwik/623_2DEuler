
//This is needed for writing/reading data from terminal
#include <iostream>
#include <string>

//This is needed to use the functions that I have defined myself in test_function_file
#include "test_function_file_header.h"


//This is necessary to define a function in the same file - the outputString function at the bottom of this file
int outputString(int num_instances);



//Well, you need a main...
int main () {

std::string test_string("This is the string to test");

//This calls the function outputString, which is defined in the same file, below.
//int something = outputString(5);

//This calls the external function return_input
//int return_input_result = return_input(50);

//This calls another external function, string_length
int string_size = string_length(test_string);


//std::cout << "return_input returns: " << return_input_result << "\n";

std::cout << "string_length returns: " << string_size << "\n";
std::cout << test_string.length() << std::endl;


return 0;
}





//This is a test function used to output a string num_instances times
int outputString(int num_instances) {

for (int i = 1; i <= num_instances; i++){
	std::cout << "Function Says Hello World!\n";
}
return (num_instances+1);
}
