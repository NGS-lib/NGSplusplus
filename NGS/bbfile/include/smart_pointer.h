/* 
 * @file   smart_pointer.h
 * @author Charles Coulombe
 *
 * @date 11 October 2011, 14:57
 *
 * Note : not thread-safe, this must be used with caution
 *
 * All class were inspired from boost::shared_ptr.
 */

#ifndef SMART_POINTER_H
#define	SMART_POINTER_H

/**
 * Delete policy for a scalar.
 * Free memory for a scalar object and set @a pointer to @c NULL.
 */
template <class TYPE>
class delete_scalar_policy
{
  public:

    /**
     * Free memory of the scalar object and set @a pointer to @c NULL.
     * @param pointer
     */
    void operator()(TYPE *&pointer) const
    {
        delete pointer;
        pointer = 0;
    }
};

/**
 * Delete policy for an array.
 * Free memory for an array and set @a pointer to @c NULL.
 */
template <class TYPE>
class delete_array_policy
{
  public:

    /**
     * Free memory of the array and set @a pointer to @c NULL.
     * @param pointer
     */
    void operator()(TYPE *&pointer) const
    {
        delete [] pointer;
        pointer = 0;
    }
};

/**
 * Reference ownership object.
 * Manage references for pointers.
 */
class reference
{
  private:
    // ----------------------- attributes
    /** reference counter */
    unsigned int _count;

  public:
    // ----------------------- methods

    /**
     * Adds a new reference.
     */
    void add()
    {
        _count++;
    }

    /**
     * Removes a reference.
     */
    void release()
    {
        _count--;
    }

    /**
     * Gets the number of active references.
     * @return number of references
     */
    unsigned int count() const
    {
        return _count;
    }
};

/**
 * Templated smart pointer implementation.
 * This implementation is quite basic but fit most needs
 * and respect the followings :
 *      - safe bool idiom;
 *      - no exception guarantees;
 *      - no cast around, make it explicit for construction
 */
template <class TYPE, template <class TYPE> class DELETER>
class smart_pointer
{
  public:
    // ----------------------- type
    /** pointer type */
    typedef TYPE pointer_type;

  private:
    // ----------------------- type
    /** destruction method */
    typedef DELETER<TYPE> delete_policy;

    /** safe bool idiom */
    typedef void (smart_pointer::*bool_type)();

    // ----------------------- attributes

    /** generic pointer */
    pointer_type * _data;

    /** reference counter */
    reference *_reference;

    // ----------------------- methods

    /**
     * Releases @c smart_pointer.
     * Memory is freed when no more active referneces exists
     */
    void _release()
    {
        _reference->release(); // release a reference

        // if no more active reference, free memory
        if(_reference->count() == 0)
        {
            delete_policy()(_data); // free data
            delete _reference; // free reference counter
        }
    }

    /**
     * Dummy function used to represent a logically @c true
     * boolean value in safe bool idiom.
     * @see http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Safe_bool
     */
    void _true() { };

  public:
    // ----------------------- constructors

    /**
     * Constructts @c smart_pointer and initializes it to @c NULL.
     */
    smart_pointer()
    {
        _data = 0;
        _reference = new reference(); // create new reference
        _reference->add(); // increment references
    }

    /**
     * Contructs @c smart_pointer from pointer to @a value
     * Increments number of references.
     * @param value source pointer
     */
    explicit smart_pointer(pointer_type *value)
    {
        _data = value; // build pointer from value
        _reference = new reference(); // create new reference
        _reference->add(); // increment references
    }

    /**
     * Copy constructor.
     * Copy data pointer and increment the reference count.
     * @param source source pointer
     */
    smart_pointer(const smart_pointer<TYPE, DELETER> &source)
    {
        _data = source._data; // copy pointer
        _reference = source._reference; // copy reference pointer
        _reference->add(); // increment references
    }

    /**
     * Destructs @c smart_pointer and free memory when there's no more references.
     * Decrements references.
     */
    ~smart_pointer()
    {
        // release a reference
        // and free memory when no more active references exists
        _release();
    }

    // ----------------------- methods

    /**
     * Dereferencing operator overload.
     * @return value
     */
    pointer_type &operator*()
    {
        return *_data;
    }

    /**
     * Indirection operator overload.
     * @return value
     */
    pointer_type *operator->()
    {
        return _data;
    }

    /**
     * Gets the raw pointer.
     * @return pointer
     */
    pointer_type *get() const
    {
        return _data;
    }

    /**
     * Swaps safely two @c smart_pointer.
     * @param other pointer to be swapped
     */
    void swap(smart_pointer<TYPE, DELETER> &other)
    {
        reference *ref = _reference;
        pointer_type *ptr = _data;

        // assign to "this"
        _data = other._data;
        _reference = other._reference;

        // assign to "other"
        other._data = ptr;
        other._reference = ref;
    }

    /**
     * Resets the current smart pointer.
     * @param pointer pointer to reset to
     */
    void reset(pointer_type *pointer = 0)
    {
        smart_pointer obj(pointer);
        swap(obj);
    }

    /**
     * Safe bool idiom.
     * Conversion to @c bool operator to ease logical tests with pointers.
     * @note returns a value that will logically be @c true if
     * <code >
     *      @c get() != 0
     * </code>
     * and else, value that is logically @c false. We don't return a real
     * @c bool to prevent unwanted implicit conversion for
     * instances where it would make no semantic sense, rather we
     * return a pointer to a member function as this will always
     * implicitly convert to @c true or @c false when used in a boolean
     * context but will not convert, for example, to an @c int type.
     * @return @c true or @c false
     * @see http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Safe_bool
     */
    operator bool_type() const
    {
        return _data ? &smart_pointer::_true : 0;
    }

    /**
     * operator= overload
     */
    smart_pointer<TYPE, DELETER> &operator=(const smart_pointer<TYPE, DELETER> &source)
    {
        // avoid self assignement
        if(this != &source)
        {
            // release a reference
            // and free memory when no more active references exists
            _release();

            // copy content
            _data = source._data; // copy pointer
            _reference = source._reference; // copy reference pointer
            _reference->add(); // increment references
        }
        return *this;
    }
};

/**
 * Helper class to ease creation of a @c smart_pointer that
 * implements a @c delete_scalar_policy .
 */
template <class TYPE>
class scalar_smart_pointer
{
  private:
    scalar_smart_pointer(); // disallow construction

  public:
    /** generic scalar type definition */
    typedef smart_pointer<TYPE, delete_scalar_policy> type;
};

/**
 * Helper class to ease creation of a @c smart_pointer that
 * implements a @c delete_array_policy .
 */
template <class TYPE>
class array_smart_pointer
{
  private:
    array_smart_pointer(); // disallow construction

  public:
    /** generic array type definition */
    typedef smart_pointer<TYPE, delete_array_policy> type;
};

#endif	/* SMART_POINTER_H */

