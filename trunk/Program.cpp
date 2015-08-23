#define VAR_P 25
#define VAR_Q 25
#define VAR_NODES 100

#include <iostream>
#include <vector>
#include <map>
#include <ctime> 
#include <cstdarg>
#include <set>
#include <fstream> 
#include <iomanip>
#include <limits>
#include <queue>
#include <list>
using namespace std;

bool rand(int p)
{
	int random_number = rand() % 100;
	return random_number < p;
}

class SomeException
{

};

template <class T>
class Structure;

template<typename T>
inline T max(T a, ...) {
	std::va_list ap;
	va_start(ap, a);
	T max = a;
 
	for (;;) {
		T arg = static_cast<T>(va_arg(ap, T));
		if (arg == 0)
			break;
		max = std::max(max, arg);
	}
 
	va_end(ap);
	return max;
}


class FileReader
{
private:

	ifstream fp_in; 

public:
	FileReader(string filepath)
	{		
		fp_in.open(filepath.c_str());    // open the streams			
	}

	ifstream& input_stream()
	{
		return this->fp_in;
	}

	~FileReader()
	{
		fp_in.close();
	}

};



template<class T>
class DataBaseMock
{
public:
	// liczba wierzcholkow
	int count;
	vector<int> people;
	vector<float> hours;
	vector<vector<int>> outs;
	int dest_task_id;
	map<int, bool> ids;

	DataBaseMock()
	{
		count = 0;
	}

	bool has_id(int index)
	{
		return ids[index] != NULL;
	}

	int size()
	{
		return people.size();
	}
	void resize(int count)
	{
		people.resize(count+1);
		hours.resize(count+1);
		outs.resize(count+1);
	}

	void init_from_file(string filepath)
	{
		FileReader fr(filepath);

		fr.input_stream() >> count;

		for (int i = 1; i <= count; ++i)
		{
			int index;
			fr.input_stream() >> index;
			resize(index);
			ids[index] = true;
			fr.input_stream() >> people[index] >> hours[index];
		}

		for (int i = 1; i <= count; ++i)
		{
			int index;
			int outs_count;

			fr.input_stream() >> index >> outs_count;

			if (outs_count == 0)
			{
				this->dest_task_id = index;
			}

			for (int j = 1; j <= outs_count; ++j)
			{
				int out_no;
				fr.input_stream() >> out_no;
				outs[index].push_back(out_no);
			}			
		}

	}

	void write_to_db(Structure<T> &st, int i = 0)
	{
		try {
			for (; i <= st.destination->get_task_number(); ++i)
			{
				if (st[i] == NULL) continue;
				
				ids[i] = true;
				count++;

				Node *node = st[i];
				int id = node->get_task_number();
				int n = node->task->get_people();
				float h = node->task->get_hours();

				if (rand(VAR_Q))
				{
					throw SomeException();
				}

				resize(id);
				this->hours[id] = h;
				this->people[id] = n;

				if (node->outs.size() == 0)
				{
					this->dest_task_id = id;
				}

				for (set<Node*, Node::compare_by_i>::iterator iter = node->outs.begin();
					iter != node->outs.end(); ++iter)
				{
					this->outs[id].push_back((*iter)->get_task_number());
				}

			}
		}
		catch (SomeException ex)
		{
			cout << " \n -----> SomeException caught (2) \n\n";
			write_to_db(st, i);
		}
	}

	void write()
	{
		cout << count << endl;
		for (int i = 0; i <= this->dest_task_id; i++)
		{
			if (ids[i] == NULL) continue;

			cout << i << " " << people[i] << " " << hours[i] << endl;
		}
		for (int i = 0; i <= this->dest_task_id; i++)
		{
			if (ids[i] == NULL) continue;

			cout << i << " " << outs[i].size();

			for (vector<int>::iterator it = outs[i].begin(); it != outs[i].end(); ++it)
			{
				cout << " " << (*it);
			}

			cout << endl;
			
		}
	}

	Structure<T>* create_structure()
	{
		Structure<TaskImpl> *st = new Structure<TaskImpl>();
		map<int,Node*> new_nodes;

		for (int i = 1; i <= size(); ++i)
		{
			try {
				if (has_id(i)) 
				{
					Node *new_node = new Node(i, people[i], hours[i]);
					new_nodes[i] = new_node;
					st->AddNewTask(new_node);
				}
			} catch (SomeException ex)
			{
				cout << "\n\n ---------> Exception caught! \n\n";
				--i;
				continue;
			}
		}
		for (int i = 1; i <= size(); ++i)
		{
			if (new_nodes[i] == NULL) continue;

			int ni = new_nodes[i]->get_task_number();
			if (outs[ni].size() > 0)
			for (int j = 0; j < outs[ni].size(); ++j)
			{
				int out_id = outs[ni][j];
				Node *out_node = new_nodes[out_id];
				st->AddTask(new_nodes[i],out_node);
			}
		}

		st->source = new_nodes[1];
		st->destination = new_nodes[dest_task_id];
		
		return st;
	}
};


class ITask {
public:
	virtual int get_task_number() = 0;
	virtual int get_people() = 0;
	virtual float get_hours() = 0;
	
};


class TaskImpl : public ITask {
private:
	int i;
	int n;
	float h;

public:
	
	TaskImpl(int _i, int _n, float _h):i(_i),n(_n),h(_h){}

	int get_task_number()
	{
		return this->i;
	}

	int get_people()
	{
		return this->n;
	}

	float get_hours()
	{
		return this->h;
	}

};




class Node
{
public:
	ITask* task;
	Node(){}
	Node(ITask* t)
	{
		this->task = t;
	}
	~Node()
	{
		delete(task);
	}
	Node(int i, int n, int h)
	{
		TaskImpl* task_impl = new TaskImpl(i, n, h);
		this->task = task_impl;
	}
	struct compare_by_i
	{
	  bool operator()(Node* s1, Node* s2) const
	  {
		  return (s1->get_task_number()) < (s2->get_task_number());
	  }
	};
	int get_task_number() {
		return this->task->get_task_number();
	}
	set<Node*,compare_by_i> outs;	
	set<Node*,compare_by_i> ins;	
};

template <class T>
class Structure
{
public:	
	Structure(){};

	Node *source;
	Node *destination;
	map<int, Node*> nodes;

	~Structure()
	{
		delete(source);
		delete(destination);
	}

	Structure(const Structure& structure)
	{
		this->source = structure.source;
		this->destination = structure.destination;
		
		map<int, Node*>::iterator iter = structure.nodes.begin();
		for (; iter != structure.nodes.end(); ++iter)
		{
			this->nodes[*(iter->first)] = *(iter->second);
		}
	}


	int size()
	{
		return this->nodes.size();
	}

	void AddTask(T *task)
	{
		Node* new_node = new Node(task);
		this->AddTask(new_node);
	}

	void AddTask(Node *node)
	{
		this->nodes[node->get_task_number()] = node;
	}

	void AddTask(Node *node, Node *node_out)
	{
		this->nodes[node->get_task_number()] = node;
		this->nodes[node_out->get_task_number()] = node_out;
		node->outs.insert(node_out);
		node_out->ins.insert(node);
	}

	void AddTask(vector<Node*> &nodes)
	{
		for (int i = 0; i < nodes.size(); ++i)
		{
			AddTask(nodes[i]);
		}
	}

	void AddNewTask(Node *node)
	{
		this->nodes[node->get_task_number()] = node;

		if (rand(VAR_P))
		{
			throw SomeException();
		}
	}

	Node* operator[](int  i) { 
		return this->nodes[i];
	}

	void apply_reduced_path(list<Node*> *path)
	{
		for (int i = 0; i < this->destination->get_task_number(); ++i)
		{
			if (nodes[i] == NULL) continue;

			if (nodes[i]->ins.size() == 0) {
				nodes[i] = NULL;
			}
		}

		for (list<Node*>::iterator iter = path->begin(); iter != path->end(); ++iter)
		{
			this->AddTask(*iter);
		}
	}
};

template <class T>
class Edge
{
public:
	T* start;
	T* end;
	Edge (T* start, T* end)
	{
		this->start = start;
		this->end = end;
	}
};


Edge<Node>* get_random_edge(vector<Edge<Node>*> *edges)
{
	int size = edges->size();

	int random_index = rand() % size;

	return (*edges)[random_index];

}

template<class T>
vector<Node*> generate_random_graph(int size)
{
	// zrodlo
	Node* source = new Node(1, rand() % 100, rand() % 100);
	// ujscie
	Node* dest = new Node(size, rand() % 100, rand() % 100);

	vector<Node*> nodes;
	nodes.push_back(source);
	nodes.push_back(dest);

	vector<Edge<Node>*> edges;

	edges.push_back(new Edge<Node>(source, dest));

	for (int i = 2; i < size; ++i)
	{
		Edge<Node>* random_edge = get_random_edge(&edges);	

		cout << "wylosowana krawedz: (" << random_edge->start->get_task_number()
			<< ", " << random_edge->end->get_task_number() << ")\n";

		Node* newide = new Node(i, rand() % 100, rand() % 100);
		
		// usuwamy stare powiazanie 
		random_edge->start->outs.erase(random_edge->end);
		random_edge->end->ins.erase(random_edge->start);

		// dodajemy nowe dwa powiazania
		random_edge->start->outs.insert(newide);
		newide->ins.insert(random_edge->start);
		
		newide->outs.insert(random_edge->end);
		random_edge->end->ins.insert(newide);

		edges.push_back(new Edge<Node>(random_edge->start, newide));
		edges.push_back(new Edge<Node>(newide, random_edge->end));

		nodes.push_back(newide);
	}

	// wypisywanie

	cout << "---- \n";

	for (int j = 0; j < nodes.size(); ++j)
	{
		set<Node*,Node::compare_by_i>::iterator iter = nodes[j]->outs.begin();

		for (iter; iter != nodes[j]->outs.end(); ++iter)
		{
			Node* current_out = *(iter);
			cout << "(" << nodes[j]->get_task_number() << ", " <<current_out->get_task_number() << ")\n";

		}
		
	}

	return nodes;

}

typedef set<Node*,Node::compare_by_i>::iterator SetIter;

// znajduje sciezke krytyczna i zwraca liste
// kolejnych wierzcholkow ze sciezki krytycznej
template<class T>
class CriticalPathFinder
{
	vector<T*> path;
	Structure<T> *structure;
	map<int, float> dist;
	list<Node*> critical_path;

public:

	CriticalPathFinder(Structure<T> *structure)
	{
		this->structure = structure;
	}

	int find_critical_path()
	{
		// inijcjalizacja
		int max_size = this->structure->size();
		for (int i = 1; i <= this->structure->size(); ++i)
		{
			Node* node = (*(this->structure))[i];
			if (node){
				int task_number = node->get_task_number();
				dist[task_number] = 0;
			}
		}
		
		queue<int> q;
		// podczas wrzucania do kolejki nalezy ustawic dist[i]
		int source_id = this->structure->source->get_task_number();
		q.push(source_id);
		dist[source_id] = this->structure->source->task->get_hours();

		while(!q.empty())
		{
			int current_id = q.front();
			q.pop();
			Node* current_node = (*(this->structure))[current_id];
			for (SetIter iter = current_node->outs.begin(); 
				iter != current_node->outs.end(); ++iter)
			{
				Node* current_out = *iter;
				dist[current_out->get_task_number()] = max(dist[current_out->get_task_number()],current_out->task->get_hours() + dist[current_node->get_task_number()]);
				q.push(current_out->get_task_number());
			}
			
		}

		cout << "path finder finished\n\n";

		// teraz mamy juz obliczone dystanse,
		// trzeba odpowiednio wybrac wierzcholki do sciezki krytycznej

		Node* current = this->structure->destination;
		while(current != NULL)
		{
	
			critical_path.push_back(current);

			if (current->ins.size() == 0) {
				current = NULL;
				continue;
			}
	
			int next_id = (*(current->ins.begin()))->get_task_number();
			for (SetIter iter = current->ins.begin(); iter != current->ins.end(); ++iter)
			{
				if (dist[(*iter)->get_task_number()] > dist[next_id])
				{
					next_id = (*iter)->get_task_number();
				}
			}

			current = (*(this->structure))[next_id];

		}

		this->critical_path.reverse();
		cout << "finished\n";

	return 0;

	}

	list<Node*> &get_critical_path()
	{
		return this->critical_path;
	}


};
 
// bierze vector z wezlami skladajacymi sie na sciezke
// i skleja pierwsze dwa wezly, ktorych sklejenie
// zmniejsza dlugosc sciezki krytycznej
class PathReducer
{
	list<Node*> path;

public:
	PathReducer(list<Node*> path)
	{
		this->path = path;
	}

	bool reduce()
	{
		list<Node*>::iterator current = this->path.begin();
			
		list<Node*>::iterator next = current;
		for(++next; next != this->path.end(); ++next) 
		{
			if (is_reduction_possible(*current, *next))
			{
				cout << "\t\tredukujemy!\n\n";
				reduce(current);
				return true;
			}
			
			current = next;
		}
		return false;
	}

	list<Node*> &get_reduced_path()
	{
		return this->path;
	}

private:
	
	bool is_reduction_possible(Node* first, Node* second)
	{
	
		if (first->ins.size() == 0 && second->outs.size() == 0) {
			return false;
		}

		int n1 = first->task->get_people();
		int n2 = second->task->get_people();
		float h1 = first->task->get_hours();
		float h2 = second->task->get_hours();
		float h = (n1*h1+n2*h2)/(n1+n2);

		return (n1+n2)* h < n1*h1 + n2*h2;
	}


	void reduce(list<Node*>::iterator current_iter)
	{
		Node* f1 = *current_iter;
		current_iter++;
		Node* f2 = *(current_iter);
		current_iter--;
		Node* new_node = new Node();
		int n1 = f1->task->get_people();
		int n2 = f2->task->get_people();
		float h1 = f1->task->get_hours();
		float h2 = f2->task->get_hours();
		float h = (n1*h1+n2*h2)/(n1+n2);
		new_node->task = new TaskImpl(f1->task->get_task_number(),n1+n2,h);

		merge(f1, f2, new_node);

		this->path.insert(current_iter, new_node);
		
		this->path.remove(f1);
		this->path.remove(f2);
		cout << "done\n";
	}

	void merge(Node *f1, Node *f2, Node *new_node)
	{
		// dla kazdego ina wchodzacego do f1
		// usuwamy outa f1 i dodajemy outa new_node
		for(set<Node*,Node::compare_by_i>::iterator iter = f1->ins.begin(); iter != f1->ins.end(); ++iter)
		{
			Node *in_node = *iter;
			in_node->outs.erase(f1);
			in_node->outs.insert(new_node);
			new_node->ins.insert(in_node);
		}

		// dla kazdego outa wychodzacego z f2
		// usuwamy ina f2 i dodajemy ina new_node
		for(set<Node*,Node::compare_by_i>::iterator iter = f2->outs.begin(); iter != f2->outs.end(); ++iter)
		{
			Node *out_node = *iter;
			out_node->ins.erase(f2);
			out_node->ins.insert(new_node);
			new_node->outs.insert(out_node);
		}

		f2->ins.erase(f1);

		// dla kazdego outa wychodzacego z f1
		// usuwamy ina z powiazanego noda i dodajemy ina new_node
		for(set<Node*,Node::compare_by_i>::iterator iter = f1->outs.begin(); iter != f1->outs.end(); ++iter)
		{
			if ((*iter)->get_task_number() == f2->get_task_number()) continue;
			Node *out_node = *iter;
			out_node->ins.erase(f1);
			out_node->ins.insert(new_node);
			new_node->outs.insert(out_node);
		}

		new_node->outs.erase(f2);
	}
};

int main()
{

	Structure<TaskImpl> random;

	srand((unsigned)time(0)); 

	vector<Node*> nodes = generate_random_graph<TaskImpl>(VAR_NODES);
	
	random.AddTask(nodes);
	random.source = nodes[0];
	random.destination = nodes[1];

/*
	st.AddTask(nodes);

	st.source = nodes[0];
	st.destination = nodes[nodes.size() - 1];

	CriticalPathFinder<TaskImpl> cpf(&st);
	cpf.find_critical_path();
	list<Node*> path = cpf.get_critical_path();

	PathReducer pr(path);

	pr.reduce();
*/

	DataBaseMock<TaskImpl> dbm;
	dbm.write_to_db(random);
//	dbm.init_from_file("c:/in3");

	bool reduced = true;
	while(reduced)
	{

		Structure<TaskImpl> *st = dbm.create_structure();

		CriticalPathFinder<TaskImpl> cpf(st);
		cpf.find_critical_path();
		list<Node*> path = cpf.get_critical_path();

		PathReducer pr(path);

		reduced = pr.reduce();

		if (!reduced) break;

		list<Node*> reduced_path = pr.get_reduced_path();

		st->apply_reduced_path(&reduced_path);

		DataBaseMock<TaskImpl> dbWriter;
		dbm = dbWriter;
		dbm.write_to_db(*st);
		dbm.write();

	}

//	test_sets();

	//system("pause");
	return 0;


}