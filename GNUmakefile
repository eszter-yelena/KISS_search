# Compiler and flags
CC = g++
#CFLAGS = -pthread -O3 -std=c++11
CFLAGS = -pthread -g -std=c++11 -O3

# Source directory and files
SOURCE_DIR = .
CPP_SOURCES = main.cpp kmerUtilities.cpp serialiseKmersMap.cpp smithWaterman.cpp hash.cpp ssw_cpp.cpp
C_SOURCES = ssw.c
SOURCES = $(CPP_SOURCES) $(C_SOURCES)

# Header directory
HEADER_DIR = headers
HEADERS = $(wildcard $(HEADER_DIR)/*.h)

# Object directory and files
OBJECT_DIR = objects
CPP_OBJECTS = $(CPP_SOURCES:%.cpp=$(OBJECT_DIR)/%.o)
C_OBJECTS = $(C_SOURCES:%.c=$(OBJECT_DIR)/%.o)
OBJECTS = $(CPP_OBJECTS) $(C_OBJECTS)

# Executable name
EXECUTABLE = KISS.out

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@

$(CPP_OBJECTS): $(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(HEADERS)
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

$(C_OBJECTS): $(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

