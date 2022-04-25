//
// Created by jins on 10/20/20.
//

#ifndef MAXINF_BRGSTKNN_2020_V2PINDEX_H
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <set>
#include <unordered_set>



class V2PNVDHashIndex{
    public:
        V2PNVDHashIndex(){}
        void setV2PIndex(unordered_map<int, int>& idx){
            this-> v2pHashMap = idx;
        }
        void setNVDIndex(unordered_map<int, int>& idx2){
            this-> nvdHashMap = idx2;
        }

        void initialIndexGiven(unordered_map<int, int>& v2p_idx, unordered_map<int, int>& nvd_idx){
            v2p_idx = v2pHashMap;
            nvd_idx = nvdHashMap;
        }
    unordered_map<int, int> v2pHashMap;
    //unordered_map<int, vector<int>> l2pHashMap;
    unordered_map<int, int> nvdHashMap;


    private:
        
        friend class boost::serialization::access;

        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->v2pHashMap;
            ar & this->nvdHashMap; // Serialized to help Dijkstra search (could easily be reconstructed from leafOccurenceList)

        }
};


class HybridPOINNHashIndex{
public:
    HybridPOINNHashIndex(){}

    void setV2PIndex(unordered_map<int, int>& idx){
        this-> v2pHashMap = idx;
    }

    void setNVDIndex(unordered_map<int, int>& idx2){
        this-> nvdHashMap = idx2;
    }

    void setL2PIndex(unordered_map<int, vector<int>>& idx3){
        this -> l2pHashMap = idx3;
    }


    void initialIndexGiven(unordered_map<int, int>& v2p_idx, unordered_map<int, vector<int>>& l2p_idx, unordered_map<int, int>& nvd_idx){
        v2p_idx =v2pHashMap;
        l2p_idx = l2pHashMap;
        nvd_idx = nvdHashMap;
    }
    unordered_map<int, int> v2pHashMap;
    unordered_map<int, vector<int>> l2pHashMap;
    unordered_map<int, int> nvdHashMap;


private:

    friend class boost::serialization::access;

    // Boost Serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->v2pHashMap;
        ar & this->l2pHashMap;
        ar & this->nvdHashMap; // Serialized to help Dijkstra search (could easily be reconstructed from leafOccurenceList)

    }
};


class SocialGraphBinary{
public:

    unordered_map<int,vector<int>> _wholeSocialGraph;
    vector<int> user_nodes;
    unordered_map<int,vector<int>> _friendshipMap;   // user,  follower
    unordered_map<int,vector<int>> _followedMap;
    //unordered_map<int,<unordered_map<int,int>>> _linkTable;

    SocialGraphBinary(){}

    void setGraphData(unordered_map<int, vector<int>>& graph_s, vector<int>& nodes_s){
        this-> _wholeSocialGraph = graph_s;
        this -> user_nodes = nodes_s;
    }

    /*void setSocialLink(unordered_map<int,vector<int>>& frindMap, unordered_map<int,vector<int>>& followMap,
                       unordered_map<int,<unordered_map<int,int>>>& linkTable){
        this->_friendshipMap = frindMap;
        this-> _followedMap = followMap;
        this-> _linkTable = linkTable;
    }*/
    void setSocialLink(unordered_map<int,vector<int>>& frindMap, unordered_map<int,vector<int>>& followMap){
        this->_friendshipMap = frindMap;
        this-> _followedMap = followMap;
    }

private:
    friend class boost::serialization::access;

    // Boost Serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->_wholeSocialGraph;
        ar & this->user_nodes;
        ar & this->_friendshipMap;
        ar & this->_followedMap;
       // ar & this->_linkTable;

    }



};

/*
class CheckInBinary{
public:

    unordered_map<int,vector<int>> _userCheckInIDList;
    unordered_map<int,vector<int>> _userCheckInCountList;

    unordered_map<int,vector<int>> _poiCheckInIDList;   // poi, checkin user...
    unordered_map<int,vector<int>> _poiCheckInCountList;

    CheckInBinary(){}

    void setUserCheckIn(unordered_map<int, vector<int>>& ucheck_pid, unordered_map<int, vector<int>>& ucheck_pcount){
        this-> _userCheckInIDList = ucheck_pid;
        this -> _userCheckInCountList = ucheck_pcount;
    }


    void setPOICheckIn(unordered_map<int,vector<int>>& pcheck_uid, unordered_map<int, vector<int>>& pcheck_ucount){
        this->_poiCheckInIDList = pcheck_uid;
        this-> _poiCheckInCountList = pcheck_ucount;
    }

private:
    friend class boost::serialization::access;

    // Boost Serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->_userCheckInIDList;
        ar & this->_userCheckInCountList;
        ar & this->_poiCheckInIDList;
        ar & this->_poiCheckInCountList;
        // ar & this->_linkTable;

    }



};
*/


class CheckInBinary{
public:

    unordered_map<int,vector<pair<int,int>>> _userCheckInfoList;  // u_id, {<poi_id, count>...}

    unordered_map<int,vector<pair<int,int>>> _poiCheckInfoList;   // poi, {<user_id, count>...}

    CheckInBinary(){}

    void setCheckInData( unordered_map<int,vector<pair<int,int>>>& ucheck,  unordered_map<int,vector<pair<int,int>>>& pcheck){
        this-> _userCheckInfoList = ucheck;
        this -> _poiCheckInfoList = pcheck;
    }



private:
    friend class boost::serialization::access;

    // Boost Serialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & this->_userCheckInfoList;
        ar & this->_poiCheckInfoList;

    }



};


#define MAXINF_BRGSTKNN_2020_V2PINDEX_H

#endif //MAXINF_BRGSTKNN_2020_V2PINDEX_H
