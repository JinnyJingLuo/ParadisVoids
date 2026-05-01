// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef CIRCULARLINKEDLIST_H_
#define CIRCULARLINKEDLIST_H_

#include "list"

using namespace std;

template <typename Type> class CircularLinkedList
{
	public:
		CircularLinkedList()
		{
			Reset();
		}
		CircularLinkedList(const CircularLinkedList& oList)
		{
			*this = oList;
		}
		~CircularLinkedList()
		{
			Reset();
		}
		CircularLinkedList& operator=(const CircularLinkedList& oList)
		{
			Reset();
			m_loList = oList.m_loList;
			ResetIterator();
			return *this;
		}
		void Reset()
		{
			m_loList.clear();
			ResetIterator();
		}
		void PreInsert(Type oItem)
		{
			m_loList.insert(m_liCurrentItem,oItem);
		}
		void PostInsert(Type oItem)
		{
			IncrementIterator();
			m_loList.insert(m_liCurrentItem,oItem);
			DecrementIterator();
		}
		void Append(Type oItem)
		{
			m_loList.push_back(oItem);
		}
		void IncrementIterator()
		{
			m_liCurrentItem++;
			if(m_liCurrentItem == m_loList.end())
			{
				m_liCurrentItem = m_loList.begin();
			}
		}
		void DecrementIterator()
		{
			if(m_liCurrentItem == m_loList.begin())
			{
				m_liCurrentItem = m_loList.end();
			}
			m_liCurrentItem--;
		}
		void ResetIterator()
		{
			m_liCurrentItem = m_loList.begin();
		}
		bool IsAtBeginning() const
		{
			return (m_liCurrentItem == m_loList.begin());
		}
		const list<Type>& GetList() const
		{
			return m_loList;
		}
		Type GetCurrentItem() const
		{
			return (*m_liCurrentItem);
		}
		Type GetPreviousItem()
		{
			DecrementIterator();
			Type oItem = (*m_liCurrentItem);
			IncrementIterator();
			return oItem;
		}
		Type GetNextItem()
		{
			IncrementIterator();
			Type oItem = (*m_liCurrentItem);
			DecrementIterator();
			return oItem;
		}
		Type* GetCurrentItemPointer()
		{
			return &(*m_liCurrentItem);
		}
		Type* GetPreviousItemPointer()
		{
			DecrementIterator();
			Type* poItem = &(*m_liCurrentItem);
			IncrementIterator();
			return poItem;
		}
		Type* GetNextItemPointer()
		{
			IncrementIterator();
			Type* poItem = &(*m_liCurrentItem);
			DecrementIterator();
			return poItem;
		}
		void DropItem()
		{
			if(!m_loList.empty())
			{
				m_liCurrentItem = m_loList.erase(m_liCurrentItem);
			}
		}
		void Reverse()
		{
			m_loList.reverse();
			ResetIterator();
		}
		unsigned int GetSize() const
		{
			return (unsigned int)m_loList.size();
		}
		bool IsEmpty() const
		{
			if(m_loList.size() == 0)
			{
				return true;
			}
			return false;
		}
	private:
	
	protected:
		list<Type> m_loList;
		typename list<Type>::iterator m_liCurrentItem;
};

#endif



