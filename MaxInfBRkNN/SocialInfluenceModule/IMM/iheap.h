
const int VectorDefaultSize = 100;
template <typename _Key, typename _Value>
struct Key_Value
{
    _Key key;
    _Value value;
    Key_Value( const _Key& k , const _Value& v )
    {
        key = k;
        value = v;
    }
    bool operator==( const Key_Value<_Key,_Value>& p )const
    {
        return ( key == p.key );
    }
    bool operator<( const Key_Value<_Key,_Value>& p )const
    {
        if ( key < p.key ) return true;
        else if ( key > p.key ) return false;
        return ( value < p.value );
    }

    Key_Value( int tmp )
    {
        key = tmp;
        value = -1;
    }

    Key_Value()
    {
    }
};

template <typename _T>
class iVector
{
    public:
        unsigned int m_size;
        _T* m_data;
        unsigned int m_num;

        void free_mem()
        {
            delete[] m_data;
        }

        iVector()
        {
            //printf("%d\n",VectorDefaultSize);
            m_size = VectorDefaultSize;
            m_data = new _T[VectorDefaultSize];
            m_num = 0;
        }
        iVector( unsigned int n )
        {
            if ( n == 0 )
            {
                n = VectorDefaultSize;
            }
            //		printf("iVector allocate: %d\n",n);
            m_size = n;
            m_data = new _T[m_size];
            m_num = 0;
        }
        void push_back( _T d )
        {
            if ( m_num == m_size )
            {
                re_allocate( m_size*2 );
            }
            m_data[m_num] = d ;
            m_num++;		
        }
        void push_back( const _T* p, unsigned int len )
        {
            while ( m_num + len > m_size )
            {
                re_allocate( m_size*2 );
            }
            memcpy( m_data+m_num, p, sizeof(_T)*len );
            m_num += len;
        }

        void re_allocate( unsigned int size )
        {
            if ( size < m_num )
            {
                return;
            }
            _T* tmp = new _T[size];
            memcpy( tmp, m_data, sizeof(_T)*m_num );
            m_size = size;
            delete[] m_data;
            m_data = tmp;
        }
        void Sort()
        {
            if ( m_num < 20 )
            {
                int k ;
                _T tmp;
                for ( int i = 0 ; i < m_num-1 ; ++i )
                {
                    k = i ;
                    for ( int j = i+1 ; j < m_num ; ++j )
                        if ( m_data[j] < m_data[k] ) k = j ;
                    if ( k != i )
                    {
                        tmp = m_data[i];
                        m_data[i] = m_data[k];
                        m_data[k] = tmp;
                    }
                }
            }
            else sort( m_data, m_data+m_num );
        }
        void unique()
        {
            if ( m_num == 0 ) return;
            Sort();
            unsigned int j = 0;
            for ( unsigned int i = 0 ; i < m_num ; ++i )
                if ( !(m_data[i] == m_data[j]) )
                {
                    ++j;
                    if ( j != i ) m_data[j] = m_data[i];
                }
            m_num = j+1;
        }
        int BinarySearch( _T& data )
        {
            for ( int x = 0 , y = m_num-1 ; x <= y ; )
            {
                int p = (x+y)/2;
                if ( m_data[p] == data ) return p;
                if ( m_data[p] < data ) x = p+1;
                else y = p-1;
            }
            return -1;
        }
        void clean()
        {
            m_num = 0;
        }
        void assign( iVector& t )
        {
            m_num = t.m_num;
            m_size = t.m_size;
            delete[] m_data;
            m_data = t.m_data;
        }

        bool remove( _T& x )
        {
            for ( int l = 0 , r = m_num ; l < r ; )
            {
                int m = (l+r)/2;

                if ( m_data[m] == x )
                {
                    m_num--;
                    if ( m_num > m ) memmove( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                    return true;
                }
                else if ( m_data[m] < x ) l = m+1;
                else r = m;
            }
            return false;
        }

        void sorted_insert( _T& x )
        {
            if ( m_num == 0 )
            {
                push_back( x );
                return;
            }

            if ( m_num == m_size ) re_allocate( m_size*2 );

            int l,r;

            for ( l = 0 , r = m_num ; l < r ; )
            {
                int m = (l+r)/2;
                if ( m_data[m] < x ) l = m+1;
                else r = m;
            }

            if ( l < m_num && m_data[l] == x )
            {
                printf("Insert Duplicate....\n");
                //		break;
            }
            else
            {
                if ( m_num > l )
                {
                    memmove( m_data+l+1, m_data+l, sizeof(_T)*(m_num-l) );
                }
                m_num++;
                m_data[l] = x;
            }
        }

        bool remove_unsorted( _T& x )
        {
            for ( int m = 0 ; m < m_num ; ++m )
            {
                if ( m_data[m] == x )
                {
                    m_num--;
                    if ( m_num > m ) memcpy( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                    return true;
                }
            }
            return false;
        }

        _T& operator[]( unsigned int i )
        {
            //if ( i < 0 || i >= m_num ) 
            //{
            //	printf("iVector [] out of range!!!\n");
            //}
            return m_data[i];
        }
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //close range check for [] in iVector if release

};
template <typename _T>
struct iMap
{	
    _T* m_data;
    int m_num;
    int cur;
    iVector<int> occur;

    iMap()
    {
        m_data = NULL;
        m_num = 0;
    }

    void free_mem()
    {
        delete[] m_data;
        occur.free_mem();
    }

    void initialize( int n )
    {
        occur.clean();
        m_num = n;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = 1073741834;
        cur = 0;
    }
    void clean()
    {
        for ( int i = 0 ; i <(int) occur.m_num ; ++i )
        {
            m_data[occur[i]] = 1073741834;
        }
        occur.clean();
        cur = 0;
    }
    _T get( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //	printf("iMap get out of range!!!\n");
        //	return -8;
        //}
        return m_data[p];
    }
    void erase( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //	printf("iMap get out of range!!!\n");
        //}
        m_data[p] = 1073741834;
        cur--;
    }
    bool notexist( int p )
    {
        return m_data[p] == 1073741834 ;
    }
    bool exist( int p )
    {
        return !(m_data[p] == 1073741834);
    }
    void insert( int p , _T d )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //	printf("iMap insert out of range!!!\n");
        //}
        if ( m_data[p] == 1073741834 )
        {
            occur.push_back( p );
            cur++;
        }
        m_data[p] = d;
    }
    void inc( int p )
    {
        //if ( m_data[p] == 1073741834 )
        //{
        //	printf("inc some unexisted point\n");
        //}
        m_data[p]++;
    }
    void inc( int p , int x )
    {
        //if ( m_data[p] == 1073741834 )
        //{
        //	printf("inc some unexisted point\n");
        //}
        m_data[p] += x;
    }
    void dec( int p )
    {
        //if ( m_data[p] == 1073741834 )
        //{
        //	printf("dec some unexisted point\n" );
        //}
        m_data[p]--;
    }
    //close range check when release!!!!!!!!!!!!!!!!!!!!	
};


template <typename _Value>
struct iHeap
{
    iMap<int> pos;
    iVector<Key_Value<int,_Value> > m_data;

    bool inited = false;
    void initialize( int n )
    {
        inited = true;
        pos.initialize(n);		
        m_data.re_allocate(1000);
    }

    int head()
    {
        return m_data[0].key;
    }

    void smart_clean(vector<int>& clean_helper){
        m_data.clean();
        for(auto element: clean_helper)
            pos.erase(element);
    }

    void clean()
    {
        pos.occur.clean();
    }

    void DeepClean()
    {
        pos.clean();
        m_data.clean();
    }

    void insert( int x, _Value y )
    {
        if ( pos.notexist(x) )
        {
            m_data.push_back( Key_Value<int,_Value>(x,y) );
            pos.insert( x , m_data.m_num-1 );
            up( m_data.m_num-1 );
        }
        else
        {
            if ( y < m_data[pos.get(x)].value )
            {
                m_data[ pos.get(x) ].value = y;
                up( pos.get(x) );
            }
            else
            {
                m_data[ pos.get(x) ].value = y;
                down( pos.get(x) );
            }
        }
    }
    int pop()
    {
        int tmp = m_data[0].key;
        pos.erase(tmp);
        if ( m_data.m_num == 1 ) 
        {
            m_data.m_num--;
            return tmp;
        }
        m_data[0] = m_data[m_data.m_num-1];
        m_data.m_num--;
        down(0);
        return tmp;
    }
    bool empty()
    {
        return ( m_data.m_num == 0 );
    }
    void up( int p )
    {
        Key_Value<int,_Value> x = m_data[p];
        for ( ; p > 0 && x.value < m_data[(p-1)/2].value ; p = (p-1)/2 )
        {
            m_data[p] = m_data[(p-1)/2];
            pos.insert(m_data[p].key, p);
        }
        m_data[p] = x;
        pos.insert(x.key,p);
    }
    void down( int p )
    {
        //Key_Value<int,int> x = m_data[m_data.m_num-1];
        //m_data.m_num--;

        Key_Value<int,_Value> tmp;

        for ( int i ; p < (int) m_data.m_num ; p = i )
        {
            //m_data[p] = x;
            if ( p*2+1 < (int)m_data.m_num && m_data[p*2+1].value < m_data[p].value )
                i = p*2+1;
            else i = p;
            if ( p*2+2 < (int)m_data.m_num && m_data[p*2+2].value < m_data[i].value )
                i = p*2+2;
            if ( i == p ) break;

            tmp = m_data[p];
            m_data[p] = m_data[i];
            pos.insert( m_data[p].key, p );
            m_data[i] = tmp;
        }

        pos.insert(m_data[p].key,p);				
    }

};
