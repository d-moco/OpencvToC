﻿# CMakeList.txt: OpencvToC 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("OpencvToC")

# 将源代码添加到此项目的可执行文件。
add_executable (OpencvToC "OpencvToC.cpp" "OpencvToC.h" "core/cvMat.h" "core/cvMat.c" "core/cvAlloc.h" "core/imgpro.h" "core/imgpro.c" "sod_img_reader.h" "sod_img_writer.h" "sod.h" "sod.c" "core/tracking.h" "core/types_c.h" "core/boundingRect.h" "core/boundingRect.c")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET OpencvToC PROPERTY CXX_STANDARD 20)
endif()

# TODO: 如有需要，请添加测试并安装目标。
