CC = g++ -std=c++20 -O3 -Wall -Wextra -o
TCC = g++ -std=c++20 -ggdb -Wall -Wextra -o
INC = src/road_network.cpp src/util.cpp

default:  index query
all: main index query topcut test
main:
	$(CC) cut src/main.cpp $(INC)
index:
	$(CC) index src/index.cpp $(INC)
query:
	$(CC) query src/query.cpp $(INC)
generate_bucket_query:
	$(CC) generate_bucket_query src/generate_queries.cpp $(INC)
idxAndQuery:
	$(CC) idxAndQuery src/idxAndQuery.cpp $(INC)
topcut:
	$(CC) topcut src/topcut.cpp $(INC)
test:
	$(TCC) test src/test.cpp $(INC)
clean:
	rm -f cut index query topcut test idxAndQuery generate_bucket_query

.PHONY: main index query topcut test idxAndQuery generate_bucket_query
 