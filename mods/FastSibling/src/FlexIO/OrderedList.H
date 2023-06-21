#ifndef __ORDEREDLIST__TMPL_
#define __ORDEREDLIST__TMPL_
/*------------------------------------------------------------
  Typical use of the orderd list type
  
   List<sphere *> SphereList;
   sphere *s;
   SphereList.reset();
   while(s=SphereList.getNext())
       s->visible=1;
------------------------------------------------------------*/
// #include <stdio.h>
// #include <iostream.h>
// #include <string.h>
// #include "NewString.H"

#define FOREACH(item,list) list.reset(); while(item=list.getNext())
#define PFOREACH(item,list) list->reset(); while(item=list->getNext())

template <class T>
class OrderedList
{
  class LLink
  {		
  public:
    T data;		
    LLink *next,*prev;
    LLink(T d,LLink *nxt,LLink *prv=0):
    data(d),next(nxt),prev(prv)
    {
      if(nxt)
	nxt->prev=this;
      if(prv)
	prv->next=this;
    }
    ~LLink()
    {
      if(prev)
	prev->next=next;
      if(next)
	next->prev=prev;	
    }
  };
  LLink *head,*current,*tail;
  int size;
public:
  T nullentry;
  OrderedList():head(0),size(0),current(0),nullentry(0){}
  OrderedList(OrderedList<T> &src):head(0),size(0),current(0),nullentry(0){
    T elem;
    FOREACH(elem,src)
      this->add(elem);
  }
  ~OrderedList(){
    while(head) {
      LLink *tmp=head->next;
      delete head;
      head=tmp;
    }
  }
  void add(T data)
  {
    head=new LLink(data,head);
    size++;
  }
  inline void addToHead(T data){add(data);}
  void addToTail(T data){
    // Whaaaaaa! Totally Bogus Man!!
    // This type aint complete!!!
  }
  void insertAfterCurrent(T data){
    if(!current || !current->next)
      addToTail(data);
    else
      current=new LLink(data,current->next,current);
  }
  void insertBeforeCurrent(T data) {
    if(!current)
      addToTail(data);
    else{
      if(!current->prev)
	add(data);
      else
	new LLink(data,current,current->prev);
    }
  }
  void remove(T data)
  {
    for(LLink *tmp=head;tmp;tmp=tmp->next) {
      if(tmp->data==data){
	if(tmp==head) head=tmp->next;
	if(tmp==current) current=tmp->next;
	delete tmp;
	size--;
      }
    }
  }
  void reset() { current=head; }
  T getNext() { 
    T rval;
    if(!current) return nullentry;
    rval=current->data;
    current=current->next;
    return rval;
  }		
  inline int getSize(){return size;}
  OrderedList &operator=(OrderedList<T> &src){
    if(this==&src) return *this;
    T elem;
    FOREACH(elem,src)
      this->add(elem);
    return *this;
  }
};

#endif // __ORDEREDLIST__TMPL_
