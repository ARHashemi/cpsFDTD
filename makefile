#
# 'make depend' uses makedepend to automatically generate dependencies
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = g++

# define any compile-time flags -g -Warray-bounds
CFLAGS = -Wall -fopenmp -std=c++17 -pedantic -Dcimg_use_vt100 -Dcimg_display=1
# -fsanitize=undefined#-fsanitize=address -fno-omit-frame-pointer# `pkg-config --cflags opencv`

#-m32 -fno-inline -fno-omit-frame-pointer
# define any directories containing header files other than /usr/include
# -I/usr/include/linux
INCLUDES = -I/media/arh/HOME/Workspace/FDTD/arhFDTD -I/media/arh/HOME/Workspace/CImg-3.1.4 -I/usr/include -I/usr/include/c++/11 -I/usr/lib/gcc/x86_64-linux-gnu/11/include/ -I/usr/include/x86_64-linux-gnu/c++/11/
#-I/home/arh/workspace/FDTD1/scr -I/usr/include/c++/7.2.0 -I/usr/include/c++/7 -I/usr/lib/gcc -I/usr/include/x86_64-linux-gnu/c++/7 -I/usr/lib/gcc/x86_64-linux-gnu/7/include -I/usr/include/x86_64-linux-gnu -I/home/arh/CImg# -I/usr/local/opencv/include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:# #-L/usr/lib  -L../lib
LFLAGS =#-fsanitize=address# `pkg-config --libs opencv`# -L/usr/lib/opencv/#

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lX11 -lpthread# -lopencv_core -lopencv_videoio -lopencv_highgui#-lmylib -lm -llibopencv_highgui.so.3.2.0

# define the C source files
SRCS = main.cpp BasicFDTD.cpp CPML.cpp Fields.cpp stdafx.cpp Structure.cpp TF_SF.cpp Ports.cpp DL_model.cpp OutImage.cpp Fluorescence.cpp

# CFS_PML.cpp Drude.cpp FourLevelLaser.cpp Gaussian.cpp SemiFluorescence.cpp Sources.cpp Drude_Lorentz.cpp
# define the C object files ThreeLevel.cpp IOimage.cpp
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file
MAIN = bin/Debug/FDTD

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
Debug:	all

.PHONY:	depend clean

all:	$(MAIN)
	@echo  Simple compiler named $(MAIN) has been compiled

$(MAIN):	$(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend:	$(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

main.o: stdafx.h /usr/include/stdio.h /usr/include/bits/libc-header-start.h
main.o: /usr/include/features.h /usr/include/features-time64.h
main.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
main.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
main.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
main.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
main.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
main.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
main.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
main.o: /usr/include/bits/types/__mbstate_t.h
main.o: /usr/include/bits/types/__fpos64_t.h /usr/include/bits/types/__FILE.h
main.o: /usr/include/bits/types/FILE.h /usr/include/bits/types/struct_FILE.h
main.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
main.o: /usr/include/bits/floatn-common.h /usr/include/sys/types.h
main.o: /usr/include/bits/types/clock_t.h /usr/include/bits/types/clockid_t.h
main.o: /usr/include/bits/types/time_t.h /usr/include/bits/types/timer_t.h
main.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
main.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
main.o: /usr/include/bits/byteswap.h /usr/include/bits/uintn-identity.h
main.o: /usr/include/sys/select.h /usr/include/bits/select.h
main.o: /usr/include/bits/types/sigset_t.h
main.o: /usr/include/bits/types/__sigset_t.h
main.o: /usr/include/bits/types/struct_timeval.h
main.o: /usr/include/bits/types/struct_timespec.h
main.o: /usr/include/bits/pthreadtypes.h
main.o: /usr/include/bits/thread-shared-types.h
main.o: /usr/include/bits/pthreadtypes-arch.h
main.o: /usr/include/bits/atomic_wide_counter.h
main.o: /usr/include/bits/struct_mutex.h /usr/include/bits/struct_rwlock.h
main.o: /usr/include/c++/11/iostream
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
main.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
main.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
main.o: /usr/include/c++/11/bits/memoryfwd.h
main.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
main.o: /usr/include/wchar.h /usr/include/bits/wchar.h
main.o: /usr/include/bits/types/wint_t.h /usr/include/bits/types/mbstate_t.h
main.o: /usr/include/bits/types/locale_t.h
main.o: /usr/include/bits/types/__locale_t.h /usr/include/c++/11/exception
main.o: /usr/include/c++/11/bits/exception.h
main.o: /usr/include/c++/11/bits/char_traits.h
main.o: /usr/include/c++/11/bits/stl_algobase.h
main.o: /usr/include/c++/11/bits/functexcept.h
main.o: /usr/include/c++/11/bits/exception_defines.h
main.o: /usr/include/c++/11/bits/cpp_type_traits.h
main.o: /usr/include/c++/11/ext/type_traits.h
main.o: /usr/include/c++/11/ext/numeric_traits.h
main.o: /usr/include/c++/11/bits/stl_pair.h /usr/include/c++/11/bits/move.h
main.o: /usr/include/c++/11/bits/concept_check.h
main.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
main.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
main.o: /usr/include/c++/11/debug/assertions.h
main.o: /usr/include/c++/11/bits/stl_iterator.h
main.o: /usr/include/c++/11/bits/ptr_traits.h
main.o: /usr/include/c++/11/debug/debug.h
main.o: /usr/include/c++/11/bits/predefined_ops.h
main.o: /usr/include/c++/11/bits/localefwd.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
main.o: /usr/include/c++/11/clocale /usr/include/locale.h
main.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
main.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
main.o: /usr/include/c++/11/ext/atomicity.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
main.o: /usr/include/pthread.h /usr/include/sched.h /usr/include/bits/sched.h
main.o: /usr/include/bits/types/struct_sched_param.h
main.o: /usr/include/bits/cpu-set.h /usr/include/time.h
main.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
main.o: /usr/include/bits/types/struct_itimerspec.h
main.o: /usr/include/bits/setjmp.h
main.o: /usr/include/bits/types/struct___jmp_buf_tag.h
main.o: /usr/include/bits/pthread_stack_min-dynamic.h
main.o: /usr/include/bits/pthread_stack_min.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
main.o: /usr/include/c++/11/bits/locale_classes.h /usr/include/c++/11/string
main.o: /usr/include/c++/11/bits/allocator.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
main.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
main.o: /usr/include/c++/11/bits/ostream_insert.h
main.o: /usr/include/c++/11/bits/cxxabi_forced.h
main.o: /usr/include/c++/11/bits/stl_function.h
main.o: /usr/include/c++/11/backward/binders.h
main.o: /usr/include/c++/11/bits/range_access.h
main.o: /usr/include/c++/11/bits/basic_string.h
main.o: /usr/include/c++/11/ext/alloc_traits.h
main.o: /usr/include/c++/11/bits/alloc_traits.h
main.o: /usr/include/c++/11/bits/stl_construct.h
main.o: /usr/include/c++/11/bits/basic_string.tcc
main.o: /usr/include/c++/11/bits/locale_classes.tcc
main.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
main.o: /usr/include/c++/11/bits/streambuf.tcc
main.o: /usr/include/c++/11/bits/basic_ios.h
main.o: /usr/include/c++/11/bits/locale_facets.h /usr/include/c++/11/cwctype
main.o: /usr/include/wctype.h /usr/include/bits/wctype-wchar.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
main.o: /usr/include/c++/11/bits/streambuf_iterator.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
main.o: /usr/include/c++/11/bits/locale_facets.tcc
main.o: /usr/include/c++/11/bits/basic_ios.tcc
main.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
main.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
main.o: /usr/include/math.h /usr/include/bits/math-vector.h
main.o: /usr/include/bits/libm-simd-decl-stubs.h
main.o: /usr/include/bits/flt-eval-method.h /usr/include/bits/fp-logb.h
main.o: /usr/include/bits/fp-fast.h
main.o: /usr/include/bits/mathcalls-helper-functions.h
main.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathcalls-narrow.h
main.o: /usr/include/bits/iscanonical.h /usr/include/c++/11/iomanip
main.o: /usr/include/c++/11/fstream /usr/include/c++/11/bits/codecvt.h
main.o: /usr/include/c++/11/cstdio
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
main.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
main.o: /usr/include/c++/11/bits/fstream.tcc /usr/include/c++/11/cerrno
main.o: /usr/include/errno.h /usr/include/bits/errno.h
main.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
main.o: /usr/include/asm-generic/errno.h
main.o: /usr/include/asm-generic/errno-base.h /usr/include/string.h
main.o: /usr/include/strings.h /usr/lib/gcc/x86_64-linux-gnu/11/include/omp.h
main.o: /usr/include/c++/11/ctime /media/arh/HOME/Workspace/CImg-3.1.4/CImg.h
main.o: /usr/include/c++/11/cstdlib /usr/include/stdlib.h
main.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
main.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
main.o: /usr/include/c++/11/bits/std_abs.h /usr/include/stdlib.h
main.o: /usr/include/c++/11/cstdarg /usr/include/c++/11/cstring
main.o: /usr/include/c++/11/cmath /usr/include/math.h
main.o: /usr/include/c++/11/cfloat
main.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h
main.o: /usr/include/c++/11/climits /usr/include/limits.h
main.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/limits.h
main.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
main.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
main.o: /usr/include/c++/11/algorithm /usr/include/c++/11/utility
main.o: /usr/include/c++/11/bits/stl_relops.h
main.o: /usr/include/c++/11/bits/stl_algo.h
main.o: /usr/include/c++/11/bits/algorithmfwd.h
main.o: /usr/include/c++/11/bits/stl_heap.h
main.o: /usr/include/c++/11/bits/stl_tempbuf.h /usr/include/sys/time.h
main.o: /usr/include/sys/stat.h /usr/include/bits/stat.h
main.o: /usr/include/bits/struct_stat.h /usr/include/unistd.h
main.o: /usr/include/bits/posix_opt.h /usr/include/bits/environments.h
main.o: /usr/include/bits/confname.h /usr/include/bits/getopt_posix.h
main.o: /usr/include/bits/getopt_core.h /usr/include/bits/unistd_ext.h
main.o: /usr/include/dirent.h /usr/include/bits/dirent.h
main.o: /usr/include/bits/dirent_ext.h /usr/include/fnmatch.h
main.o: /usr/include/X11/Xlib.h /usr/include/X11/X.h
main.o: /usr/include/X11/Xfuncproto.h /usr/include/X11/Xosdefs.h
main.o: /usr/include/X11/Xutil.h /usr/include/X11/keysym.h
main.o: /usr/include/X11/keysymdef.h Constants.h Fields.h Structure.h
main.o: BasicFDTD.h CPML.h OutImage.h Ports.h TF_SF.h DL_model.h
main.o: Fluorescence.h /usr/include/c++/11/random
main.o: /usr/include/c++/11/bits/c++0x_warning.h /usr/include/c++/11/complex
main.o: /usr/include/c++/11/sstream /usr/include/c++/11/bits/sstream.tcc
BasicFDTD.o: stdafx.h /usr/include/stdio.h
BasicFDTD.o: /usr/include/bits/libc-header-start.h /usr/include/features.h
BasicFDTD.o: /usr/include/features-time64.h /usr/include/bits/wordsize.h
BasicFDTD.o: /usr/include/bits/timesize.h /usr/include/stdc-predef.h
BasicFDTD.o: /usr/include/sys/cdefs.h /usr/include/bits/long-double.h
BasicFDTD.o: /usr/include/gnu/stubs.h
BasicFDTD.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
BasicFDTD.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
BasicFDTD.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
BasicFDTD.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
BasicFDTD.o: /usr/include/bits/types/__mbstate_t.h
BasicFDTD.o: /usr/include/bits/types/__fpos64_t.h
BasicFDTD.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
BasicFDTD.o: /usr/include/bits/types/struct_FILE.h
BasicFDTD.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
BasicFDTD.o: /usr/include/bits/floatn-common.h /usr/include/sys/types.h
BasicFDTD.o: /usr/include/bits/types/clock_t.h
BasicFDTD.o: /usr/include/bits/types/clockid_t.h
BasicFDTD.o: /usr/include/bits/types/time_t.h
BasicFDTD.o: /usr/include/bits/types/timer_t.h
BasicFDTD.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
BasicFDTD.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
BasicFDTD.o: /usr/include/bits/byteswap.h /usr/include/bits/uintn-identity.h
BasicFDTD.o: /usr/include/sys/select.h /usr/include/bits/select.h
BasicFDTD.o: /usr/include/bits/types/sigset_t.h
BasicFDTD.o: /usr/include/bits/types/__sigset_t.h
BasicFDTD.o: /usr/include/bits/types/struct_timeval.h
BasicFDTD.o: /usr/include/bits/types/struct_timespec.h
BasicFDTD.o: /usr/include/bits/pthreadtypes.h
BasicFDTD.o: /usr/include/bits/thread-shared-types.h
BasicFDTD.o: /usr/include/bits/pthreadtypes-arch.h
BasicFDTD.o: /usr/include/bits/atomic_wide_counter.h
BasicFDTD.o: /usr/include/bits/struct_mutex.h
BasicFDTD.o: /usr/include/bits/struct_rwlock.h /usr/include/c++/11/iostream
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
BasicFDTD.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
BasicFDTD.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
BasicFDTD.o: /usr/include/c++/11/bits/memoryfwd.h
BasicFDTD.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
BasicFDTD.o: /usr/include/wchar.h /usr/include/bits/wchar.h
BasicFDTD.o: /usr/include/bits/types/wint_t.h
BasicFDTD.o: /usr/include/bits/types/mbstate_t.h
BasicFDTD.o: /usr/include/bits/types/locale_t.h
BasicFDTD.o: /usr/include/bits/types/__locale_t.h
BasicFDTD.o: /usr/include/c++/11/exception
BasicFDTD.o: /usr/include/c++/11/bits/exception.h
BasicFDTD.o: /usr/include/c++/11/bits/char_traits.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_algobase.h
BasicFDTD.o: /usr/include/c++/11/bits/functexcept.h
BasicFDTD.o: /usr/include/c++/11/bits/exception_defines.h
BasicFDTD.o: /usr/include/c++/11/bits/cpp_type_traits.h
BasicFDTD.o: /usr/include/c++/11/ext/type_traits.h
BasicFDTD.o: /usr/include/c++/11/ext/numeric_traits.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_pair.h
BasicFDTD.o: /usr/include/c++/11/bits/move.h
BasicFDTD.o: /usr/include/c++/11/bits/concept_check.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
BasicFDTD.o: /usr/include/c++/11/debug/assertions.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_iterator.h
BasicFDTD.o: /usr/include/c++/11/bits/ptr_traits.h
BasicFDTD.o: /usr/include/c++/11/debug/debug.h
BasicFDTD.o: /usr/include/c++/11/bits/predefined_ops.h
BasicFDTD.o: /usr/include/c++/11/bits/localefwd.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
BasicFDTD.o: /usr/include/c++/11/clocale /usr/include/locale.h
BasicFDTD.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
BasicFDTD.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
BasicFDTD.o: /usr/include/c++/11/ext/atomicity.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
BasicFDTD.o: /usr/include/pthread.h /usr/include/sched.h
BasicFDTD.o: /usr/include/bits/sched.h
BasicFDTD.o: /usr/include/bits/types/struct_sched_param.h
BasicFDTD.o: /usr/include/bits/cpu-set.h /usr/include/time.h
BasicFDTD.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
BasicFDTD.o: /usr/include/bits/types/struct_itimerspec.h
BasicFDTD.o: /usr/include/bits/setjmp.h
BasicFDTD.o: /usr/include/bits/types/struct___jmp_buf_tag.h
BasicFDTD.o: /usr/include/bits/pthread_stack_min-dynamic.h
BasicFDTD.o: /usr/include/bits/pthread_stack_min.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
BasicFDTD.o: /usr/include/c++/11/bits/locale_classes.h
BasicFDTD.o: /usr/include/c++/11/string /usr/include/c++/11/bits/allocator.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
BasicFDTD.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
BasicFDTD.o: /usr/include/c++/11/bits/ostream_insert.h
BasicFDTD.o: /usr/include/c++/11/bits/cxxabi_forced.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_function.h
BasicFDTD.o: /usr/include/c++/11/backward/binders.h
BasicFDTD.o: /usr/include/c++/11/bits/range_access.h
BasicFDTD.o: /usr/include/c++/11/bits/basic_string.h
BasicFDTD.o: /usr/include/c++/11/ext/alloc_traits.h
BasicFDTD.o: /usr/include/c++/11/bits/alloc_traits.h
BasicFDTD.o: /usr/include/c++/11/bits/stl_construct.h
BasicFDTD.o: /usr/include/c++/11/bits/basic_string.tcc
BasicFDTD.o: /usr/include/c++/11/bits/locale_classes.tcc
BasicFDTD.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
BasicFDTD.o: /usr/include/c++/11/bits/streambuf.tcc
BasicFDTD.o: /usr/include/c++/11/bits/basic_ios.h
BasicFDTD.o: /usr/include/c++/11/bits/locale_facets.h
BasicFDTD.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
BasicFDTD.o: /usr/include/bits/wctype-wchar.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
BasicFDTD.o: /usr/include/c++/11/bits/streambuf_iterator.h
BasicFDTD.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
BasicFDTD.o: /usr/include/c++/11/bits/locale_facets.tcc
BasicFDTD.o: /usr/include/c++/11/bits/basic_ios.tcc
BasicFDTD.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
BasicFDTD.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
BasicFDTD.o: BasicFDTD.h Fields.h Structure.h
BasicFDTD.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/omp.h
CPML.o: stdafx.h /usr/include/stdio.h /usr/include/bits/libc-header-start.h
CPML.o: /usr/include/features.h /usr/include/features-time64.h
CPML.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
CPML.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
CPML.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
CPML.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
CPML.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
CPML.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
CPML.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
CPML.o: /usr/include/bits/types/__mbstate_t.h
CPML.o: /usr/include/bits/types/__fpos64_t.h /usr/include/bits/types/__FILE.h
CPML.o: /usr/include/bits/types/FILE.h /usr/include/bits/types/struct_FILE.h
CPML.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
CPML.o: /usr/include/bits/floatn-common.h /usr/include/sys/types.h
CPML.o: /usr/include/bits/types/clock_t.h /usr/include/bits/types/clockid_t.h
CPML.o: /usr/include/bits/types/time_t.h /usr/include/bits/types/timer_t.h
CPML.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
CPML.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
CPML.o: /usr/include/bits/byteswap.h /usr/include/bits/uintn-identity.h
CPML.o: /usr/include/sys/select.h /usr/include/bits/select.h
CPML.o: /usr/include/bits/types/sigset_t.h
CPML.o: /usr/include/bits/types/__sigset_t.h
CPML.o: /usr/include/bits/types/struct_timeval.h
CPML.o: /usr/include/bits/types/struct_timespec.h
CPML.o: /usr/include/bits/pthreadtypes.h
CPML.o: /usr/include/bits/thread-shared-types.h
CPML.o: /usr/include/bits/pthreadtypes-arch.h
CPML.o: /usr/include/bits/atomic_wide_counter.h
CPML.o: /usr/include/bits/struct_mutex.h /usr/include/bits/struct_rwlock.h
CPML.o: /usr/include/c++/11/iostream
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
CPML.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
CPML.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
CPML.o: /usr/include/c++/11/bits/memoryfwd.h
CPML.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
CPML.o: /usr/include/wchar.h /usr/include/bits/wchar.h
CPML.o: /usr/include/bits/types/wint_t.h /usr/include/bits/types/mbstate_t.h
CPML.o: /usr/include/bits/types/locale_t.h
CPML.o: /usr/include/bits/types/__locale_t.h /usr/include/c++/11/exception
CPML.o: /usr/include/c++/11/bits/exception.h
CPML.o: /usr/include/c++/11/bits/char_traits.h
CPML.o: /usr/include/c++/11/bits/stl_algobase.h
CPML.o: /usr/include/c++/11/bits/functexcept.h
CPML.o: /usr/include/c++/11/bits/exception_defines.h
CPML.o: /usr/include/c++/11/bits/cpp_type_traits.h
CPML.o: /usr/include/c++/11/ext/type_traits.h
CPML.o: /usr/include/c++/11/ext/numeric_traits.h
CPML.o: /usr/include/c++/11/bits/stl_pair.h /usr/include/c++/11/bits/move.h
CPML.o: /usr/include/c++/11/bits/concept_check.h
CPML.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
CPML.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
CPML.o: /usr/include/c++/11/debug/assertions.h
CPML.o: /usr/include/c++/11/bits/stl_iterator.h
CPML.o: /usr/include/c++/11/bits/ptr_traits.h
CPML.o: /usr/include/c++/11/debug/debug.h
CPML.o: /usr/include/c++/11/bits/predefined_ops.h
CPML.o: /usr/include/c++/11/bits/localefwd.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
CPML.o: /usr/include/c++/11/clocale /usr/include/locale.h
CPML.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
CPML.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
CPML.o: /usr/include/c++/11/ext/atomicity.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
CPML.o: /usr/include/pthread.h /usr/include/sched.h /usr/include/bits/sched.h
CPML.o: /usr/include/bits/types/struct_sched_param.h
CPML.o: /usr/include/bits/cpu-set.h /usr/include/time.h
CPML.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
CPML.o: /usr/include/bits/types/struct_itimerspec.h
CPML.o: /usr/include/bits/setjmp.h
CPML.o: /usr/include/bits/types/struct___jmp_buf_tag.h
CPML.o: /usr/include/bits/pthread_stack_min-dynamic.h
CPML.o: /usr/include/bits/pthread_stack_min.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
CPML.o: /usr/include/c++/11/bits/locale_classes.h /usr/include/c++/11/string
CPML.o: /usr/include/c++/11/bits/allocator.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
CPML.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
CPML.o: /usr/include/c++/11/bits/ostream_insert.h
CPML.o: /usr/include/c++/11/bits/cxxabi_forced.h
CPML.o: /usr/include/c++/11/bits/stl_function.h
CPML.o: /usr/include/c++/11/backward/binders.h
CPML.o: /usr/include/c++/11/bits/range_access.h
CPML.o: /usr/include/c++/11/bits/basic_string.h
CPML.o: /usr/include/c++/11/ext/alloc_traits.h
CPML.o: /usr/include/c++/11/bits/alloc_traits.h
CPML.o: /usr/include/c++/11/bits/stl_construct.h
CPML.o: /usr/include/c++/11/bits/basic_string.tcc
CPML.o: /usr/include/c++/11/bits/locale_classes.tcc
CPML.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
CPML.o: /usr/include/c++/11/bits/streambuf.tcc
CPML.o: /usr/include/c++/11/bits/basic_ios.h
CPML.o: /usr/include/c++/11/bits/locale_facets.h /usr/include/c++/11/cwctype
CPML.o: /usr/include/wctype.h /usr/include/bits/wctype-wchar.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
CPML.o: /usr/include/c++/11/bits/streambuf_iterator.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
CPML.o: /usr/include/c++/11/bits/locale_facets.tcc
CPML.o: /usr/include/c++/11/bits/basic_ios.tcc
CPML.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
CPML.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
CPML.o: Constants.h Structure.h CPML.h Fields.h /usr/include/c++/11/fstream
CPML.o: /usr/include/c++/11/bits/codecvt.h /usr/include/c++/11/cstdio
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
CPML.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
CPML.o: /usr/include/c++/11/bits/fstream.tcc /usr/include/c++/11/cerrno
CPML.o: /usr/include/errno.h /usr/include/bits/errno.h
CPML.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
CPML.o: /usr/include/asm-generic/errno.h
CPML.o: /usr/include/asm-generic/errno-base.h /usr/include/c++/11/sstream
CPML.o: /usr/include/c++/11/bits/sstream.tcc /usr/include/math.h
CPML.o: /usr/include/bits/math-vector.h
CPML.o: /usr/include/bits/libm-simd-decl-stubs.h
CPML.o: /usr/include/bits/flt-eval-method.h /usr/include/bits/fp-logb.h
CPML.o: /usr/include/bits/fp-fast.h
CPML.o: /usr/include/bits/mathcalls-helper-functions.h
CPML.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathcalls-narrow.h
CPML.o: /usr/include/bits/iscanonical.h
CPML.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/omp.h
Fields.o: stdafx.h /usr/include/stdio.h /usr/include/bits/libc-header-start.h
Fields.o: /usr/include/features.h /usr/include/features-time64.h
Fields.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
Fields.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
Fields.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
Fields.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
Fields.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
Fields.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
Fields.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
Fields.o: /usr/include/bits/types/__mbstate_t.h
Fields.o: /usr/include/bits/types/__fpos64_t.h
Fields.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
Fields.o: /usr/include/bits/types/struct_FILE.h /usr/include/bits/stdio_lim.h
Fields.o: /usr/include/bits/floatn.h /usr/include/bits/floatn-common.h
Fields.o: /usr/include/sys/types.h /usr/include/bits/types/clock_t.h
Fields.o: /usr/include/bits/types/clockid_t.h
Fields.o: /usr/include/bits/types/time_t.h /usr/include/bits/types/timer_t.h
Fields.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
Fields.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
Fields.o: /usr/include/bits/byteswap.h /usr/include/bits/uintn-identity.h
Fields.o: /usr/include/sys/select.h /usr/include/bits/select.h
Fields.o: /usr/include/bits/types/sigset_t.h
Fields.o: /usr/include/bits/types/__sigset_t.h
Fields.o: /usr/include/bits/types/struct_timeval.h
Fields.o: /usr/include/bits/types/struct_timespec.h
Fields.o: /usr/include/bits/pthreadtypes.h
Fields.o: /usr/include/bits/thread-shared-types.h
Fields.o: /usr/include/bits/pthreadtypes-arch.h
Fields.o: /usr/include/bits/atomic_wide_counter.h
Fields.o: /usr/include/bits/struct_mutex.h /usr/include/bits/struct_rwlock.h
Fields.o: /usr/include/c++/11/iostream
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
Fields.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
Fields.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
Fields.o: /usr/include/c++/11/bits/memoryfwd.h
Fields.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
Fields.o: /usr/include/wchar.h /usr/include/bits/wchar.h
Fields.o: /usr/include/bits/types/wint_t.h
Fields.o: /usr/include/bits/types/mbstate_t.h
Fields.o: /usr/include/bits/types/locale_t.h
Fields.o: /usr/include/bits/types/__locale_t.h /usr/include/c++/11/exception
Fields.o: /usr/include/c++/11/bits/exception.h
Fields.o: /usr/include/c++/11/bits/char_traits.h
Fields.o: /usr/include/c++/11/bits/stl_algobase.h
Fields.o: /usr/include/c++/11/bits/functexcept.h
Fields.o: /usr/include/c++/11/bits/exception_defines.h
Fields.o: /usr/include/c++/11/bits/cpp_type_traits.h
Fields.o: /usr/include/c++/11/ext/type_traits.h
Fields.o: /usr/include/c++/11/ext/numeric_traits.h
Fields.o: /usr/include/c++/11/bits/stl_pair.h /usr/include/c++/11/bits/move.h
Fields.o: /usr/include/c++/11/bits/concept_check.h
Fields.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
Fields.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
Fields.o: /usr/include/c++/11/debug/assertions.h
Fields.o: /usr/include/c++/11/bits/stl_iterator.h
Fields.o: /usr/include/c++/11/bits/ptr_traits.h
Fields.o: /usr/include/c++/11/debug/debug.h
Fields.o: /usr/include/c++/11/bits/predefined_ops.h
Fields.o: /usr/include/c++/11/bits/localefwd.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
Fields.o: /usr/include/c++/11/clocale /usr/include/locale.h
Fields.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
Fields.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
Fields.o: /usr/include/c++/11/ext/atomicity.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
Fields.o: /usr/include/pthread.h /usr/include/sched.h
Fields.o: /usr/include/bits/sched.h
Fields.o: /usr/include/bits/types/struct_sched_param.h
Fields.o: /usr/include/bits/cpu-set.h /usr/include/time.h
Fields.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
Fields.o: /usr/include/bits/types/struct_itimerspec.h
Fields.o: /usr/include/bits/setjmp.h
Fields.o: /usr/include/bits/types/struct___jmp_buf_tag.h
Fields.o: /usr/include/bits/pthread_stack_min-dynamic.h
Fields.o: /usr/include/bits/pthread_stack_min.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
Fields.o: /usr/include/c++/11/bits/locale_classes.h
Fields.o: /usr/include/c++/11/string /usr/include/c++/11/bits/allocator.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
Fields.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
Fields.o: /usr/include/c++/11/bits/ostream_insert.h
Fields.o: /usr/include/c++/11/bits/cxxabi_forced.h
Fields.o: /usr/include/c++/11/bits/stl_function.h
Fields.o: /usr/include/c++/11/backward/binders.h
Fields.o: /usr/include/c++/11/bits/range_access.h
Fields.o: /usr/include/c++/11/bits/basic_string.h
Fields.o: /usr/include/c++/11/ext/alloc_traits.h
Fields.o: /usr/include/c++/11/bits/alloc_traits.h
Fields.o: /usr/include/c++/11/bits/stl_construct.h
Fields.o: /usr/include/c++/11/bits/basic_string.tcc
Fields.o: /usr/include/c++/11/bits/locale_classes.tcc
Fields.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
Fields.o: /usr/include/c++/11/bits/streambuf.tcc
Fields.o: /usr/include/c++/11/bits/basic_ios.h
Fields.o: /usr/include/c++/11/bits/locale_facets.h
Fields.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
Fields.o: /usr/include/bits/wctype-wchar.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
Fields.o: /usr/include/c++/11/bits/streambuf_iterator.h
Fields.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
Fields.o: /usr/include/c++/11/bits/locale_facets.tcc
Fields.o: /usr/include/c++/11/bits/basic_ios.tcc
Fields.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
Fields.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
Fields.o: Fields.h Structure.h
stdafx.o: stdafx.h /usr/include/stdio.h /usr/include/bits/libc-header-start.h
stdafx.o: /usr/include/features.h /usr/include/features-time64.h
stdafx.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
stdafx.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
stdafx.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
stdafx.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
stdafx.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
stdafx.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
stdafx.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
stdafx.o: /usr/include/bits/types/__mbstate_t.h
stdafx.o: /usr/include/bits/types/__fpos64_t.h
stdafx.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
stdafx.o: /usr/include/bits/types/struct_FILE.h /usr/include/bits/stdio_lim.h
stdafx.o: /usr/include/bits/floatn.h /usr/include/bits/floatn-common.h
stdafx.o: /usr/include/sys/types.h /usr/include/bits/types/clock_t.h
stdafx.o: /usr/include/bits/types/clockid_t.h
stdafx.o: /usr/include/bits/types/time_t.h /usr/include/bits/types/timer_t.h
stdafx.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
stdafx.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
stdafx.o: /usr/include/bits/byteswap.h /usr/include/bits/uintn-identity.h
stdafx.o: /usr/include/sys/select.h /usr/include/bits/select.h
stdafx.o: /usr/include/bits/types/sigset_t.h
stdafx.o: /usr/include/bits/types/__sigset_t.h
stdafx.o: /usr/include/bits/types/struct_timeval.h
stdafx.o: /usr/include/bits/types/struct_timespec.h
stdafx.o: /usr/include/bits/pthreadtypes.h
stdafx.o: /usr/include/bits/thread-shared-types.h
stdafx.o: /usr/include/bits/pthreadtypes-arch.h
stdafx.o: /usr/include/bits/atomic_wide_counter.h
stdafx.o: /usr/include/bits/struct_mutex.h /usr/include/bits/struct_rwlock.h
stdafx.o: /usr/include/c++/11/iostream
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
stdafx.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
stdafx.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
stdafx.o: /usr/include/c++/11/bits/memoryfwd.h
stdafx.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
stdafx.o: /usr/include/wchar.h /usr/include/bits/wchar.h
stdafx.o: /usr/include/bits/types/wint_t.h
stdafx.o: /usr/include/bits/types/mbstate_t.h
stdafx.o: /usr/include/bits/types/locale_t.h
stdafx.o: /usr/include/bits/types/__locale_t.h /usr/include/c++/11/exception
stdafx.o: /usr/include/c++/11/bits/exception.h
stdafx.o: /usr/include/c++/11/bits/char_traits.h
stdafx.o: /usr/include/c++/11/bits/stl_algobase.h
stdafx.o: /usr/include/c++/11/bits/functexcept.h
stdafx.o: /usr/include/c++/11/bits/exception_defines.h
stdafx.o: /usr/include/c++/11/bits/cpp_type_traits.h
stdafx.o: /usr/include/c++/11/ext/type_traits.h
stdafx.o: /usr/include/c++/11/ext/numeric_traits.h
stdafx.o: /usr/include/c++/11/bits/stl_pair.h /usr/include/c++/11/bits/move.h
stdafx.o: /usr/include/c++/11/bits/concept_check.h
stdafx.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
stdafx.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
stdafx.o: /usr/include/c++/11/debug/assertions.h
stdafx.o: /usr/include/c++/11/bits/stl_iterator.h
stdafx.o: /usr/include/c++/11/bits/ptr_traits.h
stdafx.o: /usr/include/c++/11/debug/debug.h
stdafx.o: /usr/include/c++/11/bits/predefined_ops.h
stdafx.o: /usr/include/c++/11/bits/localefwd.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
stdafx.o: /usr/include/c++/11/clocale /usr/include/locale.h
stdafx.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
stdafx.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
stdafx.o: /usr/include/c++/11/ext/atomicity.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
stdafx.o: /usr/include/pthread.h /usr/include/sched.h
stdafx.o: /usr/include/bits/sched.h
stdafx.o: /usr/include/bits/types/struct_sched_param.h
stdafx.o: /usr/include/bits/cpu-set.h /usr/include/time.h
stdafx.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
stdafx.o: /usr/include/bits/types/struct_itimerspec.h
stdafx.o: /usr/include/bits/setjmp.h
stdafx.o: /usr/include/bits/types/struct___jmp_buf_tag.h
stdafx.o: /usr/include/bits/pthread_stack_min-dynamic.h
stdafx.o: /usr/include/bits/pthread_stack_min.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
stdafx.o: /usr/include/c++/11/bits/locale_classes.h
stdafx.o: /usr/include/c++/11/string /usr/include/c++/11/bits/allocator.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
stdafx.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
stdafx.o: /usr/include/c++/11/bits/ostream_insert.h
stdafx.o: /usr/include/c++/11/bits/cxxabi_forced.h
stdafx.o: /usr/include/c++/11/bits/stl_function.h
stdafx.o: /usr/include/c++/11/backward/binders.h
stdafx.o: /usr/include/c++/11/bits/range_access.h
stdafx.o: /usr/include/c++/11/bits/basic_string.h
stdafx.o: /usr/include/c++/11/ext/alloc_traits.h
stdafx.o: /usr/include/c++/11/bits/alloc_traits.h
stdafx.o: /usr/include/c++/11/bits/stl_construct.h
stdafx.o: /usr/include/c++/11/bits/basic_string.tcc
stdafx.o: /usr/include/c++/11/bits/locale_classes.tcc
stdafx.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
stdafx.o: /usr/include/c++/11/bits/streambuf.tcc
stdafx.o: /usr/include/c++/11/bits/basic_ios.h
stdafx.o: /usr/include/c++/11/bits/locale_facets.h
stdafx.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
stdafx.o: /usr/include/bits/wctype-wchar.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
stdafx.o: /usr/include/c++/11/bits/streambuf_iterator.h
stdafx.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
stdafx.o: /usr/include/c++/11/bits/locale_facets.tcc
stdafx.o: /usr/include/c++/11/bits/basic_ios.tcc
stdafx.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
stdafx.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
Structure.o: stdafx.h /usr/include/stdio.h
Structure.o: /usr/include/bits/libc-header-start.h /usr/include/features.h
Structure.o: /usr/include/features-time64.h /usr/include/bits/wordsize.h
Structure.o: /usr/include/bits/timesize.h /usr/include/stdc-predef.h
Structure.o: /usr/include/sys/cdefs.h /usr/include/bits/long-double.h
Structure.o: /usr/include/gnu/stubs.h
Structure.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
Structure.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
Structure.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
Structure.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
Structure.o: /usr/include/bits/types/__mbstate_t.h
Structure.o: /usr/include/bits/types/__fpos64_t.h
Structure.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
Structure.o: /usr/include/bits/types/struct_FILE.h
Structure.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
Structure.o: /usr/include/bits/floatn-common.h /usr/include/sys/types.h
Structure.o: /usr/include/bits/types/clock_t.h
Structure.o: /usr/include/bits/types/clockid_t.h
Structure.o: /usr/include/bits/types/time_t.h
Structure.o: /usr/include/bits/types/timer_t.h
Structure.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
Structure.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
Structure.o: /usr/include/bits/byteswap.h /usr/include/bits/uintn-identity.h
Structure.o: /usr/include/sys/select.h /usr/include/bits/select.h
Structure.o: /usr/include/bits/types/sigset_t.h
Structure.o: /usr/include/bits/types/__sigset_t.h
Structure.o: /usr/include/bits/types/struct_timeval.h
Structure.o: /usr/include/bits/types/struct_timespec.h
Structure.o: /usr/include/bits/pthreadtypes.h
Structure.o: /usr/include/bits/thread-shared-types.h
Structure.o: /usr/include/bits/pthreadtypes-arch.h
Structure.o: /usr/include/bits/atomic_wide_counter.h
Structure.o: /usr/include/bits/struct_mutex.h
Structure.o: /usr/include/bits/struct_rwlock.h /usr/include/c++/11/iostream
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
Structure.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
Structure.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
Structure.o: /usr/include/c++/11/bits/memoryfwd.h
Structure.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
Structure.o: /usr/include/wchar.h /usr/include/bits/wchar.h
Structure.o: /usr/include/bits/types/wint_t.h
Structure.o: /usr/include/bits/types/mbstate_t.h
Structure.o: /usr/include/bits/types/locale_t.h
Structure.o: /usr/include/bits/types/__locale_t.h
Structure.o: /usr/include/c++/11/exception
Structure.o: /usr/include/c++/11/bits/exception.h
Structure.o: /usr/include/c++/11/bits/char_traits.h
Structure.o: /usr/include/c++/11/bits/stl_algobase.h
Structure.o: /usr/include/c++/11/bits/functexcept.h
Structure.o: /usr/include/c++/11/bits/exception_defines.h
Structure.o: /usr/include/c++/11/bits/cpp_type_traits.h
Structure.o: /usr/include/c++/11/ext/type_traits.h
Structure.o: /usr/include/c++/11/ext/numeric_traits.h
Structure.o: /usr/include/c++/11/bits/stl_pair.h
Structure.o: /usr/include/c++/11/bits/move.h
Structure.o: /usr/include/c++/11/bits/concept_check.h
Structure.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
Structure.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
Structure.o: /usr/include/c++/11/debug/assertions.h
Structure.o: /usr/include/c++/11/bits/stl_iterator.h
Structure.o: /usr/include/c++/11/bits/ptr_traits.h
Structure.o: /usr/include/c++/11/debug/debug.h
Structure.o: /usr/include/c++/11/bits/predefined_ops.h
Structure.o: /usr/include/c++/11/bits/localefwd.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
Structure.o: /usr/include/c++/11/clocale /usr/include/locale.h
Structure.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
Structure.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
Structure.o: /usr/include/c++/11/ext/atomicity.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
Structure.o: /usr/include/pthread.h /usr/include/sched.h
Structure.o: /usr/include/bits/sched.h
Structure.o: /usr/include/bits/types/struct_sched_param.h
Structure.o: /usr/include/bits/cpu-set.h /usr/include/time.h
Structure.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
Structure.o: /usr/include/bits/types/struct_itimerspec.h
Structure.o: /usr/include/bits/setjmp.h
Structure.o: /usr/include/bits/types/struct___jmp_buf_tag.h
Structure.o: /usr/include/bits/pthread_stack_min-dynamic.h
Structure.o: /usr/include/bits/pthread_stack_min.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
Structure.o: /usr/include/c++/11/bits/locale_classes.h
Structure.o: /usr/include/c++/11/string /usr/include/c++/11/bits/allocator.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
Structure.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
Structure.o: /usr/include/c++/11/bits/ostream_insert.h
Structure.o: /usr/include/c++/11/bits/cxxabi_forced.h
Structure.o: /usr/include/c++/11/bits/stl_function.h
Structure.o: /usr/include/c++/11/backward/binders.h
Structure.o: /usr/include/c++/11/bits/range_access.h
Structure.o: /usr/include/c++/11/bits/basic_string.h
Structure.o: /usr/include/c++/11/ext/alloc_traits.h
Structure.o: /usr/include/c++/11/bits/alloc_traits.h
Structure.o: /usr/include/c++/11/bits/stl_construct.h
Structure.o: /usr/include/c++/11/bits/basic_string.tcc
Structure.o: /usr/include/c++/11/bits/locale_classes.tcc
Structure.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
Structure.o: /usr/include/c++/11/bits/streambuf.tcc
Structure.o: /usr/include/c++/11/bits/basic_ios.h
Structure.o: /usr/include/c++/11/bits/locale_facets.h
Structure.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
Structure.o: /usr/include/bits/wctype-wchar.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
Structure.o: /usr/include/c++/11/bits/streambuf_iterator.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
Structure.o: /usr/include/c++/11/bits/locale_facets.tcc
Structure.o: /usr/include/c++/11/bits/basic_ios.tcc
Structure.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
Structure.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
Structure.o: Structure.h Constants.h /usr/include/c++/11/fstream
Structure.o: /usr/include/c++/11/bits/codecvt.h /usr/include/c++/11/cstdio
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
Structure.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
Structure.o: /usr/include/c++/11/bits/fstream.tcc /usr/include/c++/11/cerrno
Structure.o: /usr/include/errno.h /usr/include/bits/errno.h
Structure.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
Structure.o: /usr/include/asm-generic/errno.h
Structure.o: /usr/include/asm-generic/errno-base.h
Structure.o: /usr/include/c++/11/sstream /usr/include/c++/11/bits/sstream.tcc
Structure.o: /usr/include/math.h /usr/include/bits/math-vector.h
Structure.o: /usr/include/bits/libm-simd-decl-stubs.h
Structure.o: /usr/include/bits/flt-eval-method.h /usr/include/bits/fp-logb.h
Structure.o: /usr/include/bits/fp-fast.h
Structure.o: /usr/include/bits/mathcalls-helper-functions.h
Structure.o: /usr/include/bits/mathcalls.h
Structure.o: /usr/include/bits/mathcalls-narrow.h
Structure.o: /usr/include/bits/iscanonical.h
TF_SF.o: TF_SF.h Structure.h /usr/include/c++/11/string
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
TF_SF.o: /usr/include/features.h /usr/include/features-time64.h
TF_SF.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
TF_SF.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
TF_SF.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
TF_SF.o: /usr/include/c++/11/bits/stringfwd.h
TF_SF.o: /usr/include/c++/11/bits/memoryfwd.h
TF_SF.o: /usr/include/c++/11/bits/char_traits.h
TF_SF.o: /usr/include/c++/11/bits/stl_algobase.h
TF_SF.o: /usr/include/c++/11/bits/functexcept.h
TF_SF.o: /usr/include/c++/11/bits/exception_defines.h
TF_SF.o: /usr/include/c++/11/bits/cpp_type_traits.h
TF_SF.o: /usr/include/c++/11/ext/type_traits.h
TF_SF.o: /usr/include/c++/11/ext/numeric_traits.h
TF_SF.o: /usr/include/c++/11/bits/stl_pair.h /usr/include/c++/11/bits/move.h
TF_SF.o: /usr/include/c++/11/bits/concept_check.h
TF_SF.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
TF_SF.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
TF_SF.o: /usr/include/c++/11/debug/assertions.h
TF_SF.o: /usr/include/c++/11/bits/stl_iterator.h
TF_SF.o: /usr/include/c++/11/bits/ptr_traits.h
TF_SF.o: /usr/include/c++/11/debug/debug.h
TF_SF.o: /usr/include/c++/11/bits/predefined_ops.h
TF_SF.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
TF_SF.o: /usr/include/wchar.h /usr/include/bits/libc-header-start.h
TF_SF.o: /usr/include/bits/floatn.h /usr/include/bits/floatn-common.h
TF_SF.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
TF_SF.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
TF_SF.o: /usr/include/bits/wchar.h /usr/include/bits/types/wint_t.h
TF_SF.o: /usr/include/bits/types/mbstate_t.h
TF_SF.o: /usr/include/bits/types/__mbstate_t.h
TF_SF.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
TF_SF.o: /usr/include/bits/types/locale_t.h
TF_SF.o: /usr/include/bits/types/__locale_t.h
TF_SF.o: /usr/include/c++/11/bits/allocator.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
TF_SF.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
TF_SF.o: /usr/include/c++/11/bits/exception.h
TF_SF.o: /usr/include/c++/11/bits/localefwd.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
TF_SF.o: /usr/include/c++/11/clocale /usr/include/locale.h
TF_SF.o: /usr/include/bits/locale.h /usr/include/c++/11/iosfwd
TF_SF.o: /usr/include/c++/11/cctype /usr/include/ctype.h
TF_SF.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
TF_SF.o: /usr/include/bits/time64.h /usr/include/bits/endian.h
TF_SF.o: /usr/include/bits/endianness.h
TF_SF.o: /usr/include/c++/11/bits/ostream_insert.h
TF_SF.o: /usr/include/c++/11/bits/cxxabi_forced.h
TF_SF.o: /usr/include/c++/11/bits/stl_function.h
TF_SF.o: /usr/include/c++/11/backward/binders.h
TF_SF.o: /usr/include/c++/11/bits/range_access.h
TF_SF.o: /usr/include/c++/11/bits/basic_string.h
TF_SF.o: /usr/include/c++/11/ext/atomicity.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
TF_SF.o: /usr/include/pthread.h /usr/include/sched.h
TF_SF.o: /usr/include/bits/types/time_t.h
TF_SF.o: /usr/include/bits/types/struct_timespec.h /usr/include/bits/sched.h
TF_SF.o: /usr/include/bits/types/struct_sched_param.h
TF_SF.o: /usr/include/bits/cpu-set.h /usr/include/time.h
TF_SF.o: /usr/include/bits/time.h /usr/include/bits/types/clock_t.h
TF_SF.o: /usr/include/bits/types/struct_tm.h
TF_SF.o: /usr/include/bits/types/clockid_t.h
TF_SF.o: /usr/include/bits/types/timer_t.h
TF_SF.o: /usr/include/bits/types/struct_itimerspec.h
TF_SF.o: /usr/include/bits/pthreadtypes.h
TF_SF.o: /usr/include/bits/thread-shared-types.h
TF_SF.o: /usr/include/bits/pthreadtypes-arch.h
TF_SF.o: /usr/include/bits/atomic_wide_counter.h
TF_SF.o: /usr/include/bits/struct_mutex.h /usr/include/bits/struct_rwlock.h
TF_SF.o: /usr/include/bits/setjmp.h /usr/include/bits/types/__sigset_t.h
TF_SF.o: /usr/include/bits/types/struct___jmp_buf_tag.h
TF_SF.o: /usr/include/bits/pthread_stack_min-dynamic.h
TF_SF.o: /usr/include/bits/pthread_stack_min.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
TF_SF.o: /usr/include/c++/11/ext/alloc_traits.h
TF_SF.o: /usr/include/c++/11/bits/alloc_traits.h
TF_SF.o: /usr/include/c++/11/bits/stl_construct.h
TF_SF.o: /usr/include/c++/11/bits/basic_string.tcc Fields.h stdafx.h
TF_SF.o: /usr/include/stdio.h /usr/include/bits/types/__fpos_t.h
TF_SF.o: /usr/include/bits/types/__fpos64_t.h
TF_SF.o: /usr/include/bits/types/struct_FILE.h /usr/include/bits/stdio_lim.h
TF_SF.o: /usr/include/sys/types.h /usr/include/bits/stdint-intn.h
TF_SF.o: /usr/include/endian.h /usr/include/bits/byteswap.h
TF_SF.o: /usr/include/bits/uintn-identity.h /usr/include/sys/select.h
TF_SF.o: /usr/include/bits/select.h /usr/include/bits/types/sigset_t.h
TF_SF.o: /usr/include/bits/types/struct_timeval.h
TF_SF.o: /usr/include/c++/11/iostream /usr/include/c++/11/ostream
TF_SF.o: /usr/include/c++/11/ios /usr/include/c++/11/exception
TF_SF.o: /usr/include/c++/11/bits/ios_base.h
TF_SF.o: /usr/include/c++/11/bits/locale_classes.h
TF_SF.o: /usr/include/c++/11/bits/locale_classes.tcc
TF_SF.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
TF_SF.o: /usr/include/c++/11/bits/streambuf.tcc
TF_SF.o: /usr/include/c++/11/bits/basic_ios.h
TF_SF.o: /usr/include/c++/11/bits/locale_facets.h /usr/include/c++/11/cwctype
TF_SF.o: /usr/include/wctype.h /usr/include/bits/wctype-wchar.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
TF_SF.o: /usr/include/c++/11/bits/streambuf_iterator.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
TF_SF.o: /usr/include/c++/11/bits/locale_facets.tcc
TF_SF.o: /usr/include/c++/11/bits/basic_ios.tcc
TF_SF.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
TF_SF.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
TF_SF.o: Constants.h /usr/include/c++/11/fstream
TF_SF.o: /usr/include/c++/11/bits/codecvt.h /usr/include/c++/11/cstdio
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
TF_SF.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
TF_SF.o: /usr/include/c++/11/bits/fstream.tcc /usr/include/c++/11/cerrno
TF_SF.o: /usr/include/errno.h /usr/include/bits/errno.h
TF_SF.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
TF_SF.o: /usr/include/asm-generic/errno.h
TF_SF.o: /usr/include/asm-generic/errno-base.h /usr/include/c++/11/sstream
TF_SF.o: /usr/include/c++/11/bits/sstream.tcc /usr/include/math.h
TF_SF.o: /usr/include/bits/math-vector.h
TF_SF.o: /usr/include/bits/libm-simd-decl-stubs.h
TF_SF.o: /usr/include/bits/flt-eval-method.h /usr/include/bits/fp-logb.h
TF_SF.o: /usr/include/bits/fp-fast.h
TF_SF.o: /usr/include/bits/mathcalls-helper-functions.h
TF_SF.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathcalls-narrow.h
TF_SF.o: /usr/include/bits/iscanonical.h
TF_SF.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/omp.h
Ports.o: Ports.h Fields.h Structure.h /usr/include/c++/11/string
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
Ports.o: /usr/include/features.h /usr/include/features-time64.h
Ports.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
Ports.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
Ports.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
Ports.o: /usr/include/c++/11/bits/stringfwd.h
Ports.o: /usr/include/c++/11/bits/memoryfwd.h
Ports.o: /usr/include/c++/11/bits/char_traits.h
Ports.o: /usr/include/c++/11/bits/stl_algobase.h
Ports.o: /usr/include/c++/11/bits/functexcept.h
Ports.o: /usr/include/c++/11/bits/exception_defines.h
Ports.o: /usr/include/c++/11/bits/cpp_type_traits.h
Ports.o: /usr/include/c++/11/ext/type_traits.h
Ports.o: /usr/include/c++/11/ext/numeric_traits.h
Ports.o: /usr/include/c++/11/bits/stl_pair.h /usr/include/c++/11/bits/move.h
Ports.o: /usr/include/c++/11/bits/concept_check.h
Ports.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
Ports.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
Ports.o: /usr/include/c++/11/debug/assertions.h
Ports.o: /usr/include/c++/11/bits/stl_iterator.h
Ports.o: /usr/include/c++/11/bits/ptr_traits.h
Ports.o: /usr/include/c++/11/debug/debug.h
Ports.o: /usr/include/c++/11/bits/predefined_ops.h
Ports.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
Ports.o: /usr/include/wchar.h /usr/include/bits/libc-header-start.h
Ports.o: /usr/include/bits/floatn.h /usr/include/bits/floatn-common.h
Ports.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
Ports.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
Ports.o: /usr/include/bits/wchar.h /usr/include/bits/types/wint_t.h
Ports.o: /usr/include/bits/types/mbstate_t.h
Ports.o: /usr/include/bits/types/__mbstate_t.h
Ports.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
Ports.o: /usr/include/bits/types/locale_t.h
Ports.o: /usr/include/bits/types/__locale_t.h
Ports.o: /usr/include/c++/11/bits/allocator.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
Ports.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
Ports.o: /usr/include/c++/11/bits/exception.h
Ports.o: /usr/include/c++/11/bits/localefwd.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
Ports.o: /usr/include/c++/11/clocale /usr/include/locale.h
Ports.o: /usr/include/bits/locale.h /usr/include/c++/11/iosfwd
Ports.o: /usr/include/c++/11/cctype /usr/include/ctype.h
Ports.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
Ports.o: /usr/include/bits/time64.h /usr/include/bits/endian.h
Ports.o: /usr/include/bits/endianness.h
Ports.o: /usr/include/c++/11/bits/ostream_insert.h
Ports.o: /usr/include/c++/11/bits/cxxabi_forced.h
Ports.o: /usr/include/c++/11/bits/stl_function.h
Ports.o: /usr/include/c++/11/backward/binders.h
Ports.o: /usr/include/c++/11/bits/range_access.h
Ports.o: /usr/include/c++/11/bits/basic_string.h
Ports.o: /usr/include/c++/11/ext/atomicity.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
Ports.o: /usr/include/pthread.h /usr/include/sched.h
Ports.o: /usr/include/bits/types/time_t.h
Ports.o: /usr/include/bits/types/struct_timespec.h /usr/include/bits/sched.h
Ports.o: /usr/include/bits/types/struct_sched_param.h
Ports.o: /usr/include/bits/cpu-set.h /usr/include/time.h
Ports.o: /usr/include/bits/time.h /usr/include/bits/types/clock_t.h
Ports.o: /usr/include/bits/types/struct_tm.h
Ports.o: /usr/include/bits/types/clockid_t.h
Ports.o: /usr/include/bits/types/timer_t.h
Ports.o: /usr/include/bits/types/struct_itimerspec.h
Ports.o: /usr/include/bits/pthreadtypes.h
Ports.o: /usr/include/bits/thread-shared-types.h
Ports.o: /usr/include/bits/pthreadtypes-arch.h
Ports.o: /usr/include/bits/atomic_wide_counter.h
Ports.o: /usr/include/bits/struct_mutex.h /usr/include/bits/struct_rwlock.h
Ports.o: /usr/include/bits/setjmp.h /usr/include/bits/types/__sigset_t.h
Ports.o: /usr/include/bits/types/struct___jmp_buf_tag.h
Ports.o: /usr/include/bits/pthread_stack_min-dynamic.h
Ports.o: /usr/include/bits/pthread_stack_min.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
Ports.o: /usr/include/c++/11/ext/alloc_traits.h
Ports.o: /usr/include/c++/11/bits/alloc_traits.h
Ports.o: /usr/include/c++/11/bits/stl_construct.h
Ports.o: /usr/include/c++/11/bits/basic_string.tcc /usr/include/string.h
Ports.o: /usr/include/strings.h /usr/include/c++/11/iostream
Ports.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
Ports.o: /usr/include/c++/11/exception /usr/include/c++/11/bits/ios_base.h
Ports.o: /usr/include/c++/11/bits/locale_classes.h
Ports.o: /usr/include/c++/11/bits/locale_classes.tcc
Ports.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
Ports.o: /usr/include/c++/11/bits/streambuf.tcc
Ports.o: /usr/include/c++/11/bits/basic_ios.h
Ports.o: /usr/include/c++/11/bits/locale_facets.h /usr/include/c++/11/cwctype
Ports.o: /usr/include/wctype.h /usr/include/bits/wctype-wchar.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
Ports.o: /usr/include/c++/11/bits/streambuf_iterator.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
Ports.o: /usr/include/c++/11/bits/locale_facets.tcc
Ports.o: /usr/include/c++/11/bits/basic_ios.tcc
Ports.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
Ports.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/fstream
Ports.o: /usr/include/c++/11/bits/codecvt.h /usr/include/c++/11/cstdio
Ports.o: /usr/include/stdio.h /usr/include/bits/types/__fpos_t.h
Ports.o: /usr/include/bits/types/__fpos64_t.h
Ports.o: /usr/include/bits/types/struct_FILE.h /usr/include/bits/stdio_lim.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
Ports.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
Ports.o: /usr/include/c++/11/bits/fstream.tcc /usr/include/c++/11/cerrno
Ports.o: /usr/include/errno.h /usr/include/bits/errno.h
Ports.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
Ports.o: /usr/include/asm-generic/errno.h
Ports.o: /usr/include/asm-generic/errno-base.h stdafx.h
Ports.o: /usr/include/sys/types.h /usr/include/bits/stdint-intn.h
Ports.o: /usr/include/endian.h /usr/include/bits/byteswap.h
Ports.o: /usr/include/bits/uintn-identity.h /usr/include/sys/select.h
Ports.o: /usr/include/bits/select.h /usr/include/bits/types/sigset_t.h
Ports.o: /usr/include/bits/types/struct_timeval.h /usr/include/c++/11/limits
Ports.o: Constants.h /usr/include/math.h /usr/include/bits/math-vector.h
Ports.o: /usr/include/bits/libm-simd-decl-stubs.h
Ports.o: /usr/include/bits/flt-eval-method.h /usr/include/bits/fp-logb.h
Ports.o: /usr/include/bits/fp-fast.h
Ports.o: /usr/include/bits/mathcalls-helper-functions.h
Ports.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathcalls-narrow.h
Ports.o: /usr/include/bits/iscanonical.h /usr/include/c++/11/sstream
Ports.o: /usr/include/c++/11/bits/sstream.tcc
DL_model.o: stdafx.h /usr/include/stdio.h
DL_model.o: /usr/include/bits/libc-header-start.h /usr/include/features.h
DL_model.o: /usr/include/features-time64.h /usr/include/bits/wordsize.h
DL_model.o: /usr/include/bits/timesize.h /usr/include/stdc-predef.h
DL_model.o: /usr/include/sys/cdefs.h /usr/include/bits/long-double.h
DL_model.o: /usr/include/gnu/stubs.h
DL_model.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
DL_model.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
DL_model.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
DL_model.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
DL_model.o: /usr/include/bits/types/__mbstate_t.h
DL_model.o: /usr/include/bits/types/__fpos64_t.h
DL_model.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
DL_model.o: /usr/include/bits/types/struct_FILE.h
DL_model.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
DL_model.o: /usr/include/bits/floatn-common.h /usr/include/sys/types.h
DL_model.o: /usr/include/bits/types/clock_t.h
DL_model.o: /usr/include/bits/types/clockid_t.h
DL_model.o: /usr/include/bits/types/time_t.h
DL_model.o: /usr/include/bits/types/timer_t.h /usr/include/bits/stdint-intn.h
DL_model.o: /usr/include/endian.h /usr/include/bits/endian.h
DL_model.o: /usr/include/bits/endianness.h /usr/include/bits/byteswap.h
DL_model.o: /usr/include/bits/uintn-identity.h /usr/include/sys/select.h
DL_model.o: /usr/include/bits/select.h /usr/include/bits/types/sigset_t.h
DL_model.o: /usr/include/bits/types/__sigset_t.h
DL_model.o: /usr/include/bits/types/struct_timeval.h
DL_model.o: /usr/include/bits/types/struct_timespec.h
DL_model.o: /usr/include/bits/pthreadtypes.h
DL_model.o: /usr/include/bits/thread-shared-types.h
DL_model.o: /usr/include/bits/pthreadtypes-arch.h
DL_model.o: /usr/include/bits/atomic_wide_counter.h
DL_model.o: /usr/include/bits/struct_mutex.h
DL_model.o: /usr/include/bits/struct_rwlock.h /usr/include/c++/11/iostream
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
DL_model.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
DL_model.o: /usr/include/c++/11/iosfwd /usr/include/c++/11/bits/stringfwd.h
DL_model.o: /usr/include/c++/11/bits/memoryfwd.h
DL_model.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
DL_model.o: /usr/include/wchar.h /usr/include/bits/wchar.h
DL_model.o: /usr/include/bits/types/wint_t.h
DL_model.o: /usr/include/bits/types/mbstate_t.h
DL_model.o: /usr/include/bits/types/locale_t.h
DL_model.o: /usr/include/bits/types/__locale_t.h
DL_model.o: /usr/include/c++/11/exception
DL_model.o: /usr/include/c++/11/bits/exception.h
DL_model.o: /usr/include/c++/11/bits/char_traits.h
DL_model.o: /usr/include/c++/11/bits/stl_algobase.h
DL_model.o: /usr/include/c++/11/bits/functexcept.h
DL_model.o: /usr/include/c++/11/bits/exception_defines.h
DL_model.o: /usr/include/c++/11/bits/cpp_type_traits.h
DL_model.o: /usr/include/c++/11/ext/type_traits.h
DL_model.o: /usr/include/c++/11/ext/numeric_traits.h
DL_model.o: /usr/include/c++/11/bits/stl_pair.h
DL_model.o: /usr/include/c++/11/bits/move.h
DL_model.o: /usr/include/c++/11/bits/concept_check.h
DL_model.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
DL_model.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
DL_model.o: /usr/include/c++/11/debug/assertions.h
DL_model.o: /usr/include/c++/11/bits/stl_iterator.h
DL_model.o: /usr/include/c++/11/bits/ptr_traits.h
DL_model.o: /usr/include/c++/11/debug/debug.h
DL_model.o: /usr/include/c++/11/bits/predefined_ops.h
DL_model.o: /usr/include/c++/11/bits/localefwd.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
DL_model.o: /usr/include/c++/11/clocale /usr/include/locale.h
DL_model.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
DL_model.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
DL_model.o: /usr/include/c++/11/ext/atomicity.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
DL_model.o: /usr/include/pthread.h /usr/include/sched.h
DL_model.o: /usr/include/bits/sched.h
DL_model.o: /usr/include/bits/types/struct_sched_param.h
DL_model.o: /usr/include/bits/cpu-set.h /usr/include/time.h
DL_model.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
DL_model.o: /usr/include/bits/types/struct_itimerspec.h
DL_model.o: /usr/include/bits/setjmp.h
DL_model.o: /usr/include/bits/types/struct___jmp_buf_tag.h
DL_model.o: /usr/include/bits/pthread_stack_min-dynamic.h
DL_model.o: /usr/include/bits/pthread_stack_min.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
DL_model.o: /usr/include/c++/11/bits/locale_classes.h
DL_model.o: /usr/include/c++/11/string /usr/include/c++/11/bits/allocator.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
DL_model.o: /usr/include/c++/11/ext/new_allocator.h /usr/include/c++/11/new
DL_model.o: /usr/include/c++/11/bits/ostream_insert.h
DL_model.o: /usr/include/c++/11/bits/cxxabi_forced.h
DL_model.o: /usr/include/c++/11/bits/stl_function.h
DL_model.o: /usr/include/c++/11/backward/binders.h
DL_model.o: /usr/include/c++/11/bits/range_access.h
DL_model.o: /usr/include/c++/11/bits/basic_string.h
DL_model.o: /usr/include/c++/11/ext/alloc_traits.h
DL_model.o: /usr/include/c++/11/bits/alloc_traits.h
DL_model.o: /usr/include/c++/11/bits/stl_construct.h
DL_model.o: /usr/include/c++/11/bits/basic_string.tcc
DL_model.o: /usr/include/c++/11/bits/locale_classes.tcc
DL_model.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
DL_model.o: /usr/include/c++/11/bits/streambuf.tcc
DL_model.o: /usr/include/c++/11/bits/basic_ios.h
DL_model.o: /usr/include/c++/11/bits/locale_facets.h
DL_model.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
DL_model.o: /usr/include/bits/wctype-wchar.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
DL_model.o: /usr/include/c++/11/bits/streambuf_iterator.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
DL_model.o: /usr/include/c++/11/bits/locale_facets.tcc
DL_model.o: /usr/include/c++/11/bits/basic_ios.tcc
DL_model.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
DL_model.o: /usr/include/c++/11/bits/istream.tcc /usr/include/c++/11/limits
DL_model.o: DL_model.h Structure.h Fields.h Constants.h /usr/include/math.h
DL_model.o: /usr/include/bits/math-vector.h
DL_model.o: /usr/include/bits/libm-simd-decl-stubs.h
DL_model.o: /usr/include/bits/flt-eval-method.h /usr/include/bits/fp-logb.h
DL_model.o: /usr/include/bits/fp-fast.h
DL_model.o: /usr/include/bits/mathcalls-helper-functions.h
DL_model.o: /usr/include/bits/mathcalls.h
DL_model.o: /usr/include/bits/mathcalls-narrow.h
DL_model.o: /usr/include/bits/iscanonical.h /usr/include/c++/11/fstream
DL_model.o: /usr/include/c++/11/bits/codecvt.h /usr/include/c++/11/cstdio
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
DL_model.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
DL_model.o: /usr/include/c++/11/bits/fstream.tcc /usr/include/c++/11/cerrno
DL_model.o: /usr/include/errno.h /usr/include/bits/errno.h
DL_model.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
DL_model.o: /usr/include/asm-generic/errno.h
DL_model.o: /usr/include/asm-generic/errno-base.h /usr/include/c++/11/sstream
DL_model.o: /usr/include/c++/11/bits/sstream.tcc
DL_model.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/omp.h
OutImage.o: OutImage.h /media/arh/HOME/Workspace/CImg-3.1.4/CImg.h
OutImage.o: /usr/include/c++/11/cstdio
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
OutImage.o: /usr/include/features.h /usr/include/features-time64.h
OutImage.o: /usr/include/bits/wordsize.h /usr/include/bits/timesize.h
OutImage.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
OutImage.o: /usr/include/bits/long-double.h /usr/include/gnu/stubs.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
OutImage.o: /usr/include/stdio.h /usr/include/bits/libc-header-start.h
OutImage.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
OutImage.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
OutImage.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
OutImage.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
OutImage.o: /usr/include/bits/types/__mbstate_t.h
OutImage.o: /usr/include/bits/types/__fpos64_t.h
OutImage.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
OutImage.o: /usr/include/bits/types/struct_FILE.h
OutImage.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
OutImage.o: /usr/include/bits/floatn-common.h /usr/include/c++/11/cstdlib
OutImage.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
OutImage.o: /usr/include/bits/waitstatus.h /usr/include/sys/types.h
OutImage.o: /usr/include/bits/types/clock_t.h
OutImage.o: /usr/include/bits/types/clockid_t.h
OutImage.o: /usr/include/bits/types/time_t.h
OutImage.o: /usr/include/bits/types/timer_t.h /usr/include/bits/stdint-intn.h
OutImage.o: /usr/include/endian.h /usr/include/bits/endian.h
OutImage.o: /usr/include/bits/endianness.h /usr/include/bits/byteswap.h
OutImage.o: /usr/include/bits/uintn-identity.h /usr/include/sys/select.h
OutImage.o: /usr/include/bits/select.h /usr/include/bits/types/sigset_t.h
OutImage.o: /usr/include/bits/types/__sigset_t.h
OutImage.o: /usr/include/bits/types/struct_timeval.h
OutImage.o: /usr/include/bits/types/struct_timespec.h
OutImage.o: /usr/include/bits/pthreadtypes.h
OutImage.o: /usr/include/bits/thread-shared-types.h
OutImage.o: /usr/include/bits/pthreadtypes-arch.h
OutImage.o: /usr/include/bits/atomic_wide_counter.h
OutImage.o: /usr/include/bits/struct_mutex.h
OutImage.o: /usr/include/bits/struct_rwlock.h /usr/include/alloca.h
OutImage.o: /usr/include/bits/stdlib-float.h
OutImage.o: /usr/include/c++/11/bits/std_abs.h /usr/include/stdlib.h
OutImage.o: /usr/include/c++/11/cstdarg /usr/include/c++/11/cstring
OutImage.o: /usr/include/string.h /usr/include/bits/types/locale_t.h
OutImage.o: /usr/include/bits/types/__locale_t.h /usr/include/strings.h
OutImage.o: /usr/include/c++/11/cmath
OutImage.o: /usr/include/c++/11/bits/cpp_type_traits.h
OutImage.o: /usr/include/c++/11/ext/type_traits.h /usr/include/math.h
OutImage.o: /usr/include/c++/11/cfloat
OutImage.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h
OutImage.o: /usr/include/c++/11/climits /usr/include/limits.h
OutImage.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/limits.h
OutImage.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
OutImage.o: /usr/include/linux/limits.h
OutImage.o: /usr/include/bits/pthread_stack_min-dynamic.h
OutImage.o: /usr/include/bits/pthread_stack_min.h
OutImage.o: /usr/include/bits/posix2_lim.h /usr/include/c++/11/ctime
OutImage.o: /usr/include/time.h /usr/include/bits/time.h
OutImage.o: /usr/include/bits/types/struct_tm.h
OutImage.o: /usr/include/bits/types/struct_itimerspec.h
OutImage.o: /usr/include/c++/11/exception
OutImage.o: /usr/include/c++/11/bits/exception.h
OutImage.o: /usr/include/c++/11/algorithm /usr/include/c++/11/utility
OutImage.o: /usr/include/c++/11/bits/stl_relops.h
OutImage.o: /usr/include/c++/11/bits/stl_pair.h
OutImage.o: /usr/include/c++/11/bits/move.h
OutImage.o: /usr/include/c++/11/bits/concept_check.h
OutImage.o: /usr/include/c++/11/bits/stl_algobase.h
OutImage.o: /usr/include/c++/11/bits/functexcept.h
OutImage.o: /usr/include/c++/11/bits/exception_defines.h
OutImage.o: /usr/include/c++/11/ext/numeric_traits.h
OutImage.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
OutImage.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
OutImage.o: /usr/include/c++/11/debug/assertions.h
OutImage.o: /usr/include/c++/11/bits/stl_iterator.h
OutImage.o: /usr/include/c++/11/bits/ptr_traits.h
OutImage.o: /usr/include/c++/11/debug/debug.h
OutImage.o: /usr/include/c++/11/bits/predefined_ops.h
OutImage.o: /usr/include/c++/11/bits/stl_algo.h
OutImage.o: /usr/include/c++/11/bits/algorithmfwd.h
OutImage.o: /usr/include/c++/11/bits/stl_heap.h
OutImage.o: /usr/include/c++/11/bits/stl_tempbuf.h
OutImage.o: /usr/include/c++/11/bits/stl_construct.h /usr/include/c++/11/new
OutImage.o: /usr/include/sys/time.h /usr/include/sys/stat.h
OutImage.o: /usr/include/bits/stat.h /usr/include/bits/struct_stat.h
OutImage.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
OutImage.o: /usr/include/bits/environments.h /usr/include/bits/confname.h
OutImage.o: /usr/include/bits/getopt_posix.h /usr/include/bits/getopt_core.h
OutImage.o: /usr/include/bits/unistd_ext.h /usr/include/dirent.h
OutImage.o: /usr/include/bits/dirent.h /usr/include/bits/dirent_ext.h
OutImage.o: /usr/include/fnmatch.h /usr/include/X11/Xlib.h
OutImage.o: /usr/include/X11/X.h /usr/include/X11/Xfuncproto.h
OutImage.o: /usr/include/X11/Xosdefs.h /usr/include/X11/Xutil.h
OutImage.o: /usr/include/X11/keysym.h /usr/include/X11/keysymdef.h
OutImage.o: /usr/include/pthread.h /usr/include/sched.h
OutImage.o: /usr/include/bits/sched.h
OutImage.o: /usr/include/bits/types/struct_sched_param.h
OutImage.o: /usr/include/bits/cpu-set.h /usr/include/bits/setjmp.h
OutImage.o: /usr/include/bits/types/struct___jmp_buf_tag.h
OutImage.o: /usr/include/c++/11/iostream /usr/include/c++/11/ostream
OutImage.o: /usr/include/c++/11/ios /usr/include/c++/11/iosfwd
OutImage.o: /usr/include/c++/11/bits/stringfwd.h
OutImage.o: /usr/include/c++/11/bits/memoryfwd.h
OutImage.o: /usr/include/c++/11/bits/postypes.h /usr/include/c++/11/cwchar
OutImage.o: /usr/include/wchar.h /usr/include/bits/wchar.h
OutImage.o: /usr/include/bits/types/wint_t.h
OutImage.o: /usr/include/bits/types/mbstate_t.h
OutImage.o: /usr/include/c++/11/bits/char_traits.h
OutImage.o: /usr/include/c++/11/bits/localefwd.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
OutImage.o: /usr/include/c++/11/clocale /usr/include/locale.h
OutImage.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
OutImage.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
OutImage.o: /usr/include/c++/11/ext/atomicity.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
OutImage.o: /usr/include/c++/11/bits/locale_classes.h
OutImage.o: /usr/include/c++/11/string /usr/include/c++/11/bits/allocator.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
OutImage.o: /usr/include/c++/11/ext/new_allocator.h
OutImage.o: /usr/include/c++/11/bits/ostream_insert.h
OutImage.o: /usr/include/c++/11/bits/cxxabi_forced.h
OutImage.o: /usr/include/c++/11/bits/stl_function.h
OutImage.o: /usr/include/c++/11/backward/binders.h
OutImage.o: /usr/include/c++/11/bits/range_access.h
OutImage.o: /usr/include/c++/11/bits/basic_string.h
OutImage.o: /usr/include/c++/11/ext/alloc_traits.h
OutImage.o: /usr/include/c++/11/bits/alloc_traits.h
OutImage.o: /usr/include/c++/11/bits/basic_string.tcc
OutImage.o: /usr/include/c++/11/bits/locale_classes.tcc
OutImage.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
OutImage.o: /usr/include/c++/11/bits/streambuf.tcc
OutImage.o: /usr/include/c++/11/bits/basic_ios.h
OutImage.o: /usr/include/c++/11/bits/locale_facets.h
OutImage.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
OutImage.o: /usr/include/bits/wctype-wchar.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
OutImage.o: /usr/include/c++/11/bits/streambuf_iterator.h
OutImage.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
OutImage.o: /usr/include/c++/11/bits/locale_facets.tcc
OutImage.o: /usr/include/c++/11/bits/basic_ios.tcc
OutImage.o: /usr/include/c++/11/bits/ostream.tcc /usr/include/c++/11/istream
OutImage.o: /usr/include/c++/11/bits/istream.tcc
Fluorescence.o: stdafx.h /usr/include/stdio.h
Fluorescence.o: /usr/include/bits/libc-header-start.h /usr/include/features.h
Fluorescence.o: /usr/include/features-time64.h /usr/include/bits/wordsize.h
Fluorescence.o: /usr/include/bits/timesize.h /usr/include/stdc-predef.h
Fluorescence.o: /usr/include/sys/cdefs.h /usr/include/bits/long-double.h
Fluorescence.o: /usr/include/gnu/stubs.h
Fluorescence.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h
Fluorescence.o: /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h
Fluorescence.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
Fluorescence.o: /usr/include/bits/time64.h /usr/include/bits/types/__fpos_t.h
Fluorescence.o: /usr/include/bits/types/__mbstate_t.h
Fluorescence.o: /usr/include/bits/types/__fpos64_t.h
Fluorescence.o: /usr/include/bits/types/__FILE.h
Fluorescence.o: /usr/include/bits/types/FILE.h
Fluorescence.o: /usr/include/bits/types/struct_FILE.h
Fluorescence.o: /usr/include/bits/stdio_lim.h /usr/include/bits/floatn.h
Fluorescence.o: /usr/include/bits/floatn-common.h /usr/include/sys/types.h
Fluorescence.o: /usr/include/bits/types/clock_t.h
Fluorescence.o: /usr/include/bits/types/clockid_t.h
Fluorescence.o: /usr/include/bits/types/time_t.h
Fluorescence.o: /usr/include/bits/types/timer_t.h
Fluorescence.o: /usr/include/bits/stdint-intn.h /usr/include/endian.h
Fluorescence.o: /usr/include/bits/endian.h /usr/include/bits/endianness.h
Fluorescence.o: /usr/include/bits/byteswap.h
Fluorescence.o: /usr/include/bits/uintn-identity.h /usr/include/sys/select.h
Fluorescence.o: /usr/include/bits/select.h /usr/include/bits/types/sigset_t.h
Fluorescence.o: /usr/include/bits/types/__sigset_t.h
Fluorescence.o: /usr/include/bits/types/struct_timeval.h
Fluorescence.o: /usr/include/bits/types/struct_timespec.h
Fluorescence.o: /usr/include/bits/pthreadtypes.h
Fluorescence.o: /usr/include/bits/thread-shared-types.h
Fluorescence.o: /usr/include/bits/pthreadtypes-arch.h
Fluorescence.o: /usr/include/bits/atomic_wide_counter.h
Fluorescence.o: /usr/include/bits/struct_mutex.h
Fluorescence.o: /usr/include/bits/struct_rwlock.h
Fluorescence.o: /usr/include/c++/11/iostream
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++config.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/os_defines.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/cpu_defines.h
Fluorescence.o: /usr/include/c++/11/ostream /usr/include/c++/11/ios
Fluorescence.o: /usr/include/c++/11/iosfwd
Fluorescence.o: /usr/include/c++/11/bits/stringfwd.h
Fluorescence.o: /usr/include/c++/11/bits/memoryfwd.h
Fluorescence.o: /usr/include/c++/11/bits/postypes.h
Fluorescence.o: /usr/include/c++/11/cwchar /usr/include/wchar.h
Fluorescence.o: /usr/include/bits/wchar.h /usr/include/bits/types/wint_t.h
Fluorescence.o: /usr/include/bits/types/mbstate_t.h
Fluorescence.o: /usr/include/bits/types/locale_t.h
Fluorescence.o: /usr/include/bits/types/__locale_t.h
Fluorescence.o: /usr/include/c++/11/exception
Fluorescence.o: /usr/include/c++/11/bits/exception.h
Fluorescence.o: /usr/include/c++/11/bits/char_traits.h
Fluorescence.o: /usr/include/c++/11/bits/stl_algobase.h
Fluorescence.o: /usr/include/c++/11/bits/functexcept.h
Fluorescence.o: /usr/include/c++/11/bits/exception_defines.h
Fluorescence.o: /usr/include/c++/11/bits/cpp_type_traits.h
Fluorescence.o: /usr/include/c++/11/ext/type_traits.h
Fluorescence.o: /usr/include/c++/11/ext/numeric_traits.h
Fluorescence.o: /usr/include/c++/11/bits/stl_pair.h
Fluorescence.o: /usr/include/c++/11/bits/move.h
Fluorescence.o: /usr/include/c++/11/bits/concept_check.h
Fluorescence.o: /usr/include/c++/11/bits/stl_iterator_base_types.h
Fluorescence.o: /usr/include/c++/11/bits/stl_iterator_base_funcs.h
Fluorescence.o: /usr/include/c++/11/debug/assertions.h
Fluorescence.o: /usr/include/c++/11/bits/stl_iterator.h
Fluorescence.o: /usr/include/c++/11/bits/ptr_traits.h
Fluorescence.o: /usr/include/c++/11/debug/debug.h
Fluorescence.o: /usr/include/c++/11/bits/predefined_ops.h
Fluorescence.o: /usr/include/c++/11/bits/localefwd.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++locale.h
Fluorescence.o: /usr/include/c++/11/clocale /usr/include/locale.h
Fluorescence.o: /usr/include/bits/locale.h /usr/include/c++/11/cctype
Fluorescence.o: /usr/include/ctype.h /usr/include/c++/11/bits/ios_base.h
Fluorescence.o: /usr/include/c++/11/ext/atomicity.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/gthr-default.h
Fluorescence.o: /usr/include/pthread.h /usr/include/sched.h
Fluorescence.o: /usr/include/bits/sched.h
Fluorescence.o: /usr/include/bits/types/struct_sched_param.h
Fluorescence.o: /usr/include/bits/cpu-set.h /usr/include/time.h
Fluorescence.o: /usr/include/bits/time.h /usr/include/bits/types/struct_tm.h
Fluorescence.o: /usr/include/bits/types/struct_itimerspec.h
Fluorescence.o: /usr/include/bits/setjmp.h
Fluorescence.o: /usr/include/bits/types/struct___jmp_buf_tag.h
Fluorescence.o: /usr/include/bits/pthread_stack_min-dynamic.h
Fluorescence.o: /usr/include/bits/pthread_stack_min.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/atomic_word.h
Fluorescence.o: /usr/include/c++/11/bits/locale_classes.h
Fluorescence.o: /usr/include/c++/11/string
Fluorescence.o: /usr/include/c++/11/bits/allocator.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++allocator.h
Fluorescence.o: /usr/include/c++/11/ext/new_allocator.h
Fluorescence.o: /usr/include/c++/11/new
Fluorescence.o: /usr/include/c++/11/bits/ostream_insert.h
Fluorescence.o: /usr/include/c++/11/bits/cxxabi_forced.h
Fluorescence.o: /usr/include/c++/11/bits/stl_function.h
Fluorescence.o: /usr/include/c++/11/backward/binders.h
Fluorescence.o: /usr/include/c++/11/bits/range_access.h
Fluorescence.o: /usr/include/c++/11/bits/basic_string.h
Fluorescence.o: /usr/include/c++/11/ext/alloc_traits.h
Fluorescence.o: /usr/include/c++/11/bits/alloc_traits.h
Fluorescence.o: /usr/include/c++/11/bits/stl_construct.h
Fluorescence.o: /usr/include/c++/11/bits/basic_string.tcc
Fluorescence.o: /usr/include/c++/11/bits/locale_classes.tcc
Fluorescence.o: /usr/include/c++/11/stdexcept /usr/include/c++/11/streambuf
Fluorescence.o: /usr/include/c++/11/bits/streambuf.tcc
Fluorescence.o: /usr/include/c++/11/bits/basic_ios.h
Fluorescence.o: /usr/include/c++/11/bits/locale_facets.h
Fluorescence.o: /usr/include/c++/11/cwctype /usr/include/wctype.h
Fluorescence.o: /usr/include/bits/wctype-wchar.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_base.h
Fluorescence.o: /usr/include/c++/11/bits/streambuf_iterator.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/ctype_inline.h
Fluorescence.o: /usr/include/c++/11/bits/locale_facets.tcc
Fluorescence.o: /usr/include/c++/11/bits/basic_ios.tcc
Fluorescence.o: /usr/include/c++/11/bits/ostream.tcc
Fluorescence.o: /usr/include/c++/11/istream
Fluorescence.o: /usr/include/c++/11/bits/istream.tcc
Fluorescence.o: /usr/include/c++/11/limits Fluorescence.h Structure.h
Fluorescence.o: Fields.h /usr/include/c++/11/random
Fluorescence.o: /usr/include/c++/11/bits/c++0x_warning.h
Fluorescence.o: /usr/include/c++/11/complex /usr/include/c++/11/cmath
Fluorescence.o: /usr/include/math.h /usr/include/c++/11/bits/std_abs.h
Fluorescence.o: /usr/include/stdlib.h /usr/include/c++/11/sstream
Fluorescence.o: /usr/include/c++/11/bits/sstream.tcc Constants.h
Fluorescence.o: /usr/include/c++/11/fstream
Fluorescence.o: /usr/include/c++/11/bits/codecvt.h /usr/include/c++/11/cstdio
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/basic_file.h
Fluorescence.o: /usr/include/x86_64-linux-gnu/c++/11/bits/c++io.h
Fluorescence.o: /usr/include/c++/11/bits/fstream.tcc
Fluorescence.o: /usr/include/c++/11/cerrno /usr/include/errno.h
Fluorescence.o: /usr/include/bits/errno.h /usr/include/linux/errno.h
Fluorescence.o: /usr/include/asm/errno.h /usr/include/asm-generic/errno.h
Fluorescence.o: /usr/include/asm-generic/errno-base.h /usr/include/stdlib.h
Fluorescence.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
Fluorescence.o: /usr/include/alloca.h /usr/include/bits/stdlib-float.h
Fluorescence.o: /usr/include/math.h /usr/include/bits/math-vector.h
Fluorescence.o: /usr/include/bits/libm-simd-decl-stubs.h
Fluorescence.o: /usr/include/bits/flt-eval-method.h
Fluorescence.o: /usr/include/bits/fp-logb.h /usr/include/bits/fp-fast.h
Fluorescence.o: /usr/include/bits/mathcalls-helper-functions.h
Fluorescence.o: /usr/include/bits/mathcalls.h
Fluorescence.o: /usr/include/bits/mathcalls-narrow.h
Fluorescence.o: /usr/include/bits/iscanonical.h
