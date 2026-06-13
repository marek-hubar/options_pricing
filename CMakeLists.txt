# 1. Require a minimum CMake version
cmake_minimum_required(VERSION 3.15)

# 2. Name your project
project(OptionPricer)

# 3. Force C++20 globally for this project
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional but recommended: Disable compiler-specific extensions 
# (forces strict standard C++ so it compiles on any machine)
set(CMAKE_CXX_EXTENSIONS OFF)

# 4. Generate the JSON file that your Language Server needs!
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 5. Tell CMake which source files to compile into the final executable
add_executable(pricer pricing.cpp)
