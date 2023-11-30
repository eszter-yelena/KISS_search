SOURCES = main.obj kmerUtilities.obj serialiseKmersMap.obj smithWaterman.obj hash.obj ssw_cpp.obj ssw.obj

CFLAGS = /I headers /EHsc /Zi
CPPFLAGS = $(CFLAGS) /std:c++20

KISS.exe : $(SOURCES)
	$(CPP) $(CPPFLAGS) /FeKISS.exe $(SOURCES)

clean:
	del *.obj KISS.exe

