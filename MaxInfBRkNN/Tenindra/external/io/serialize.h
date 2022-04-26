/* Original Work Copyright (C) 2005, 2006, 2007, 2008 Robert Geisberger, Dominik Schultes, Peter Sanders, Universitaet Karlsruhe (TH)
 * Modified Work Copyright (C) 2012 Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 * 
 * The following file is a derivative work of the code from 
 * http://sourceforge.net/projects/ntu-sp-exp/ developed for the paper below.
 * The authors have requested users of this code to cite the paper.
 * 
 * Lingkun Wu, Xiaokui Xiao, Dingxiong Deng, Gao Cong, Andy Diwen Zhu, Shuigeng Zhou
 * Shortest Path and Distance Queries on Road Networks: An Experimental Evaluation
 * PVLDB 5(5): 406-417 (2012)
 *
 * That work is in turn a derivative work of the code from
 * Contraction Hierarchies (http://algo2.iti.kit.edu/english/routeplanning.php),
 * licensed under AGPLv3. Thus this work is also licensed under that license.
 * 
 * Road Network kNN Experimental Evaluation is distributed in the hope 
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the 
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public 
 * License along with Road Network kNN Experimental Evaluation; see 
 * LICENSE.txt; if not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SERIALIZE_H
#define SERIALIZE_H

//#include <ext/hash_map>

/** Writes a given primitive value (e.g. int, double) to the given stream. */
template < typename value_type >
inline void writePrimitive(ostream& out, const value_type& v) {
	out.write((char*)&v, sizeof(v));
}

/** Reads a primitive value from the given stream. */
template < typename value_type >
inline void readPrimitive(istream& in, value_type& v) {
	in.read((char*)&v, sizeof(v));
}

/** Serializer for primitive values (e.g. int, double). */
template < typename value_type >
class PrimitiveSerializer
{
public:
	static void serialize(ostream& out, const value_type& v) { writePrimitive(out, v); }
	static void deserialize(istream& in, value_type& v) { readPrimitive(in, v); }
};

/** Serializer for complex objects that provide (de)serialize-methods. */
template < typename value_type >
class ComplexSerializer
{
public:
	static void serialize(ostream& out, const value_type& v) { v.serialize(out); }
	static void deserialize(istream& in, value_type& v) { v.deserialize(in); }
};

/** Serializer for stl::pairs. */
template < typename first_type, typename second_type,
typename firstSerializer = PrimitiveSerializer<first_type>,
typename secondSerializer = PrimitiveSerializer<second_type> >
class PairSerializer
{
private:
	typedef pair<first_type, second_type> value_type;

public:
	static void serialize(ostream& out, const value_type& v) {
		firstSerializer::serialize(out, v.first);
		secondSerializer::serialize(out, v.second);
	}

	static void deserialize(istream& in, value_type& v) {
		firstSerializer::deserialize(in, v.first);
		secondSerializer::deserialize(in, v.second);
	}
};

/** Serializer for arrays. */
template < typename value_type, typename size_type,
typename value_serializer = PrimitiveSerializer<value_type> >
class ArraySerializer
{
public:
	static void serialize(ostream& out, const value_type* a, const size_type n) {
		writePrimitive(out, n);
		for (size_type i = 0; i < n; i++) value_serializer::serialize(out, a[i]);
	}

	static value_type* deserialize(istream& in, size_type& n) {
		readPrimitive(in, n);
		value_type* a = new value_type[n];
		for (size_type i = 0; i < n; i++) value_serializer::deserialize(in, a[i]);
		return a;
	}
};

/** Serializer for stl::vectors. */
template < typename value_type, typename size_type,
typename value_serializer = PrimitiveSerializer<value_type> >
class VectorSerializer
{
public:
	static void serialize(ostream& out, const vector<value_type>& vec) {
		size_type n = (size_type)vec.size();
		writePrimitive(out, n);
		for (size_type i = 0; i < n; i++) 
		{
			//cout << i << endl;
			value_serializer::serialize(out, vec[i]);
		}
	}

	static void deserialize(istream& in, vector<value_type>& vec) {
		size_type n = 0;
		readPrimitive(in, n);
		vec.resize(n);
		for (size_type i = 0; i < n; i++) {
			//cout << i << endl;		
			value_serializer::deserialize(in, vec[i]);
		}
	}
};


/** Serializer for stl::set. */
template < typename value_type, typename size_type, 
typename value_serializer = PrimitiveSerializer<value_type>, typename compare = less<value_type>>
class SetSerializer
{
public:
	static void serialize(ostream& out, const set<value_type, compare>& s) {
		size_type n = (size_type)s.size();
		//cout << n << endl;
		writePrimitive(out, n);
		typename set <value_type, compare>::const_iterator sit = s.begin();
		for (size_type i = 0; i < n; i++) 
		{
			value_serializer::serialize(out, *sit);
			sit++;
		}
	}

	static void deserialize(istream& in, set<value_type, compare>& s) {
		size_type n = 0;
		readPrimitive(in, n);
		//cout << n << endl;
		value_type element;
		size_type i;
		for (i = 0; i < n; i++) {		
			value_serializer::deserialize(in, element);
			s.insert(element);
		}
		//cout << i << endl;
	}
};
/** Serializer for HashMaps. */
//template < typename key_type, typename data_type, typename size_type,
//           typename data_serializer = PrimitiveSerializer<data_type> >
//class HashMapSerializer
//{
//private:
//    typedef pair<key_type, data_type> value_type;
//    typedef PairSerializer< key_type, data_type,
//                            PrimitiveSerializer<key_type>,
//                            data_serializer > value_serializer;
//    typedef __gnu_cxx::hash_map<key_type, data_type> HashMap;
//    
//public:
//    static void serialize(ostream& out, const HashMap& hm) {
//        size_type n = (size_type)hm.size();
//        writePrimitive(out, n);
//        for (typename HashMap::const_iterator it = hm.begin(); it != hm.end(); it++)
//            value_serializer::serialize(out, (*it));
//        }
//
//    static void deserialize(istream& in, HashMap& hm) {
//        size_type n = 0;
//        readPrimitive(in, n);
//        value_type v;
//        for (size_type i = 0; i < n; i++) {
//            value_serializer::deserialize(in, v);
//            hm.insert(v);
//        }
//    }
//};



#endif // SERIALIZE_H
