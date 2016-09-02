//
// Created by icme-huang on 8/29/16.
//

#ifndef INC_2DFLUIDSTRUCTURE_NONCOPYABLE_H
#define INC_2DFLUIDSTRUCTURE_NONCOPYABLE_H


class NonCopyable {
protected:
    NonCopyable(){}
    ~NonCopyable(){}
private:
    NonCopyable(const NonCopyable&);
    void operator=(const NonCopyable&);
};




#endif //INC_2DFLUIDSTRUCTURE_NONCOPYABLE_H
