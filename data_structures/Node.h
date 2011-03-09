
#ifndef __NODE_H__
#define __NODE_H__

////////////////////////////////////////////////////////////////////////////////
// File    : Node.h                                                             
// Author  : Jonathan Eastep   email: jonathan.eastep@gmail.com                 
// Written : 16 February 2011
//                                                                              
// Copyright (C) 2011 Jonathan Eastep
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or            
// (at your option) any later version.                                          
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU             
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software Foundation       
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA                  
////////////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include "FCBase.h"

template <typename T>
class PtrNode {
protected:
        FCIntPtr  key;
        T*        value;

public:
        PtrNode(const FCIntPtr key, T* const value): key(key), value(value) {}

        ~PtrNode() {};

        virtual bool operator<(const PtrNode<T>& rhs) const {
                return key < rhs.key;
        }

        FCIntPtr getkey() {
                return key;
        }

        T* getvalue() {
                return value;
        }
};

class FCIntPtrNode : public PtrNode<FCIntPtr> {
public:
        FCIntPtrNode(const FCIntPtr keyandval): PtrNode<FCIntPtr>(keyandval, &key) {
                //there are some restrictions on the integer range
                //this may need to change if additional ops beyond add/remove are supported; base it on allowable user-space pointers
                //FIXME: we need asserts here!!!
                //assert( keyandval > FCBase<FCIntPtr>::_NULL_VALUE );  //if not, will screw up distinguishing ops
                //assert( -keyandval != FCBase<FCIntPtr>::_DEQ_VALUE ); //if not, will screw up distinguishing ops
                //assert( keyandval < FCBase<FCIntPtr>::_MAX_INT );     //if not, may or may not screw up skiplist sentinal nodes
        }

        ~FCIntPtrNode() {}

        FCIntPtr getkey() {
                return key;
        }

        FCIntPtr getvalue() {
                return *value;
        }
};


#endif
