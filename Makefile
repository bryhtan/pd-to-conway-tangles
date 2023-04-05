TARGET := pdToConwayTangles

HEADERS := $(shell find $(SRC_DIRS) -name *.h)
SRCS := $(shell find $(SRC_DIRS) -name *.c )
OBJS := $(addsuffix .o, $(basename $(SRCS)))
CFLAGS :=
DEBUG_FLAGS :=-DDEBUG

all: $(TARGET)

debug: CFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

%.o: %.c $(HEADERS)
	$(CC) -c $(CDEFINES) $(CFLAGS) $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS)
