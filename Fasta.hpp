/**
 * Copyright ©2022. Brent Weichel. All Rights Reserved.
 * Permission to use, copy, modify, and/or distribute this software, in whole
 * or part by any means, without express prior written agreement is prohibited.
 */
#pragma once

#include <algorithm>
#include <fstream>
#include <map>
#include <ostream>
#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

/**
 * Notes:
 *   - Requires C++14 minimum
 *
 * Todo:
 *   - Add support for bare sequence files.
 *   - Set IUB/IUPC amino/nucleid acid codes.
 *   - Add detection if neither mode is set for the acid codes.
 *   - Fill in 'Any' between end of sequence and offset for mutable access.
 *   - Add `operator[]( string )` to FastaFile to access the vector of FastaSequence
 *   - If the identifier changes for a sequence in the file, then we need that change reflected.
 */

/**
 * This class contains a single sequence and its associated identifier.
 * When assigning the identifier, control characters will be removed.
 * The sequence will be normalized by removing control characters and any
 * non-(amino/nucleic) characters.
 */
class FastaSequence
{
private:
	std::string mIdentifier;  // The sequence identifier.
	std::string mSequence;    // Sequence.

	void _copyAssign(
		const FastaSequence& other )
	{
		mIdentifier = other.mIdentifier;
		mSequence = other.mSequence;
	}

	void _moveAssign(
		FastaSequence&& other )
	{
		mIdentifier = std::move( other.mIdentifier );
		mSequence = std::move( other.mSequence );
	}

	// Remove any control characters
	void _normalizeIdentifier()
	{
		// Remove all non-print characters
		mIdentifier.erase(
			std::remove_if(
				mIdentifier.begin(), mIdentifier.end(),
				[ & ]( const unsigned char character )
				{
					return not std::isprint( character );
				} ),
			mIdentifier.end() );

		// Remove the leading '>' if present
		mIdentifier = mIdentifier.substr( '>' == mIdentifier[ 0 ] );

		// Remove any leading whitespace
		std::regex_replace( mIdentifier, std::regex( "^\\s+" ), std::string( "" ) );

		// Remove any trailing whitespace
		std::regex_replace( mIdentifier, std::regex( "\\s+$" ), std::string( "" ) );
	}

	// Returns true for valid sequence characters.
	bool _isValidSequenceCharacter(
		unsigned char character )
	{
		return std::isalpha( character ) or ( '-' == character ) or ( '*' == character );
	}

	// Remove control characters and
	// any non-(amino/nucleic) characters.
	void _normalizeSequence()
	{
		mSequence.erase(
			std::remove_if(
				mSequence.begin(), mSequence.end(),
				[ & ]( const unsigned char character )
				{
					return not _isValidSequenceCharacter( character );
				} ),
			mSequence.end() );
	}

public:
	/**
	 * Default constructor.
	 * @param identifier The header/sequence identifier to associate with the sequence. [default: ""]
	 * @param sequence The genomic sequence. [default: ""]
	 */
	FastaSequence(
		const std::string& identifier = std::string(),
		const std::string& sequence = std::string() )
	{
		mIdentifier = identifier;
		mSequence = sequence;

		_normalizeIdentifier();
		_normalizeSequence();
	}

	/**
	 * Move constructor.
	 * @param other R-Value to the FastaSequence to move to this instance.
	 */
	FastaSequence(
		FastaSequence&& other )
	{
		_moveAssign( std::move( other ) );
	}

	/**
	 * Copy constructor.
	 * @param other Const reference to the FastaSequence to copy to this instance.
	 */
	FastaSequence(
		const FastaSequence& other )
	{
		_copyAssign( other );
	}

	/**
	 * Append a character to the sequence count times.
	 * If the character to append is not a valid sequence character,
	 * then nothing will be appended to the sequence.
	 * @param character The character to append.
	 * @param count The number of times to append the character. [default: 1]
	 */
	void append(
		char character,
		size_t count = 1 )
	{
		if ( ( 0 < count ) and _isValidSequenceCharacter( character ) )
		{
			mSequence.append( count, character );
		}
	}

	/**
	 * Append characters from a c-string up to count characters
	 * or up to the null character of the string, which ever comes first.
	 * @param sequence The c-string to append characters from.
	 * @param count The max number of characters to append. [default: -1]
	 */
	void append(
		const char* sequence,
		size_t count = -1 )
	{
		if ( ( 0 < count ) and ( nullptr != sequence ) )
		{
			size_t length = strlen( sequence );
			count = ( count < length ) ? count : length;
			mSequence.append( sequence, count );
		}
	}

	/**
	 * Append characters from a std::string up to count characters.
	 * @param sequence The std::string to append characters from.
	 * @param count The max number of characters to append. [default: -1]
	 */
	void append(
		const std::string& sequence,
		size_t count = -1 )
	{
		this->append( sequence.c_str(), count );
	}

	/**
	 * Get the identifier for the sequence.
	 * @return A const reference to the sequence identifier.
	 */
	const std::string& identifier() const
	{
		return mIdentifier;
	}

	/**
	 * Get the length of the sequence.
	 * @return The length of the sequence is returned.
	 */
	size_t length() const
	{
		return mSequence.length();
	}

	/**
	 * The lexicographical ordering of sequences is first by the identifier, then the sequence.
	 * @param rhs The right-hand side of the comparison.
	 * @return True is returned if this instance goes before the RHS.
	 */
	bool operator<(
		const FastaSequence& rhs ) const noexcept
	{
		return ( mIdentifier == rhs.mIdentifier )
			? mSequence < rhs.mSequence
			: mIdentifier < rhs.mIdentifier;
	}

	/**
	 * The lexicographical ordering of sequences is first by the identifier, then the sequence.
	 * @param rhs The right-hand side of the comparison.
	 * @return True is returned if this instance goes before, or is equal to the RHS.
	 */
	bool operator<=(
		const FastaSequence& rhs ) const noexcept
	{
		return ( mIdentifier == rhs.mIdentifier )
			? mSequence <= rhs.mSequence
			: mIdentifier < rhs.mIdentifier;
	}

	/**
	 * The lexicographical ordering of sequences is first by the identifier, then the sequence.
	 * @param rhs The right-hand side of the comparison.
	 * @return True is returned if this instance goes after the RHS.
	 */
	bool operator>(
		const FastaSequence& rhs ) const noexcept
	{
		return ( mIdentifier == rhs.mIdentifier )
			? mSequence > rhs.mSequence
			: mIdentifier > rhs.mIdentifier;
	}

	/**
	 * The lexicographical ordering of sequences is first by the identifier, then the sequence.
	 * @param rhs The right-hand side of the comparison.
	 * @return True is returned if this instance goes after, or is equal to the RHS.
	 */
	bool operator>=(
		const FastaSequence& rhs ) const noexcept
	{
		return ( mIdentifier == rhs.mIdentifier )
			? mSequence >= rhs.mSequence
			: mIdentifier > rhs.mIdentifier;
	}

	/**
	 * Check for equality of both the identifier
	 * and the sequence; returning true if both are equal.
	 * @param other Const reference to the FastaSequence to compare against for equality.
	 * @return True is returned if both the identifier and the sequence are equal.
	 */
	bool operator==(
		const FastaSequence& other ) const
	{
		return ( mIdentifier == other.mIdentifier ) and ( mSequence == other.mSequence );
	}

	/**
	 * Check for inequality of either the identifier
	 * or the sequence; returning true if either are not equal.
	 * @param other Const reference to the FastaSequence to compare against for inequality.
	 * @return True is returned if either the identifier or the sequence are not equal.
	 */
	bool operator!=(
		const FastaSequence& other ) const
	{
		return not this->operator==( other );
	}

	/**
	 * Assign both the identifier and the sequence.
	 * @param identifierAndSequence Const pair for both the identifier and the sequence.
	 * @return Reference to this FastaSequence instance.
	 */
	FastaSequence& operator=(
		const std::pair< std::string, std::string >& identifierAndSequence )
	{
		mIdentifier = identifierAndSequence.first;
		mSequence = identifierAndSequence.second;

		_normalizeIdentifier();
		_normalizeSequence();
	}

	/**
	 * Copy assignment operator.
	 * @param other Const reference to the FastaSequence to copy to this instance.
	 * @return Reference to this FastaSequence instance.
	 */
	FastaSequence& operator=(
		const FastaSequence& other )
	{
		if ( this != &other )
		{
			_copyAssign( other );
		}

		return *this;
	}

	/**
	 * Move assignment operator.
	 * @param other R-Value to the FastaSequence to move to this instance.
	 * @return Reference to this FastaSequence instance.
	 */
	FastaSequence& operator=(
		FastaSequence&& other )
	{
		if ( this != &other )
		{
			_moveAssign( std::move( other ) );
		}

		return *this;
	}

	/**
	 * Sequence string element access.
	 * @param index Offset from the beginning of the sequence string.
	 * @return Reference to the character at the given offset.
	 */
	char& operator[](
		size_t index )
	{
		return mSequence[ index ];
	}

	/**
	 * Const sequence string element access.
	 * @param index Offset from the beginning of the sequence string.
	 * @return Const reference to the character at the given offset.
	 * @throw std::out_of_range is thrown if the index is greater than the sequence length.
	 */
	const char& operator[](
		size_t index ) const
	{
		return mSequence.at( index );
	}

	/**
	 * Get the sequence.
	 * @return A const reference to the sequence.
	 */
	const std::string& sequence() const
	{
		return mSequence;
	}

	/**
	 * Set the identifier and normalize.
	 * @param identifier Const reference to the identifier to assign to this object.
	 */
	void setIdentifier(
		const std::string& identifier )
	{
		mIdentifier = identifier;
		_normalizeIdentifier();
	}

	/**
	 * Set the sequence and normalize.
	 * @param sequence Const reference to the sequence to assign to this object.
	 */
	void setSequence(
		const std::string& sequence )
	{
		mSequence = sequence;
		_normalizeSequence();
	}
};

/**
 * This class contains a collection of FastaSequences and
 * facilitates the reading/writing of FastA files, searchinng
 * based on the identifier, iterating over the collection, addition
 * and subtraction of sequences
 */
class FastaFile
{
private:
	using IdentifierSequenceMapType = std::map< std::string, std::vector< FastaSequence > >;

	bool mIsBareSequence;
	bool mDuplicateIdentifiersAllowed;
	std::set< std::string > mIdentifiersSet;
	IdentifierSequenceMapType mIdentifierSequenceMap;

	void _copyAssign(
		const FastaFile& other )
	{
		mIsBareSequence = other.mIsBareSequence;
		mDuplicateIdentifiersAllowed = other.mDuplicateIdentifiersAllowed;
		mIdentifiersSet = other.mIdentifiersSet;
		mIdentifierSequenceMap = other.mIdentifierSequenceMap;
	}

	void _moveAssign(
		FastaFile&& other )
	{
		mIsBareSequence = std::exchange( other.mIsBareSequence, false );
		mDuplicateIdentifiersAllowed = std::exchange( other.mDuplicateIdentifiersAllowed, true );
		mIdentifiersSet = std::move( other.mIdentifiersSet );
		mIdentifierSequenceMap = std::move( other.mIdentifierSequenceMap );
	}

public:
	/**
	 * Class for iterating over the container and possibly mutating elements.
	 */
	class iterator
	{
	private:
		friend class FastaFile;

		bool mIsValid;
		IdentifierSequenceMapType::iterator mMapIterator;
		size_t mVectorOffset;

		iterator(
			IdentifierSequenceMapType::iterator mapIterator )
		{
			mIsValid = true;
			mMapIterator = mapIterator;
			mVectorOffset = 0;
		}

		void _copyAssign(
			const iterator& other )
		{
			mIsValid = other.mIsValid;
			mMapIterator = other.mMapIterator;
			mVectorOffset = other.mVectorOffset;
		}

		void _moveAssign(
			iterator&& other )
		{
			mIsValid = std::exchange( other.mIsValid, false );
			mMapIterator = std::move( other.mMapIterator );
			mVectorOffset = std::exchange( other.mVectorOffset, 0 );
		}

	public:
		using iterator_category = std::forward_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = FastaSequence;
		using pointer           = FastaSequence*;
		using reference         = FastaSequence&;

		/**
		 * Default constructor to an empty iterator.
		 */
		iterator()
		{
			mIsValid = false;
			mVectorOffset = 0;
		}

		/**
		 * Copy constructor.
		 * @param other Const reference to iterator to copy to this instance.
		 */
		iterator(
			const iterator& other )
		{
			_copyAssign( other );
		}

		/**
		 * Move constructor.
		 * @param other R-Value to iterator to move to this instance.
		 */
		iterator(
			iterator&& other )
		{
			_moveAssign( std::move( other ) );
		}

		/**
		 * Copy assignment operator.
		 * @param other Const reference to iterator to copy to this instance.
		 * @return Reference to this iterator instance is returned.
		 */
		iterator& operator=(
			const iterator& other )
		{
			if ( this != &other )
			{
				_copyAssign( other );
			}

			return *this;
		}

		/**
		 * Move assignment operator.
		 * @param other R-Value to iterator to move to this instance.
		 * @return Reference to this iterator instance is returned.
		 */
		iterator operator=(
			iterator&& other )
		{
			if ( this != &other )
			{
				_moveAssign( std::move( other ) );
			}

			return *this;
		}

		/**
		 * Compare iterators for equality.
		 * @param other Const reference to the iterator to compare against for equality.
		 * @return True is returned if this and other are equal.
		 */
		bool operator==(
			const iterator& other )
		{
			return mIsValid
				and ( mIsValid == other.mIsValid )
				and ( mMapIterator == other.mMapIterator )
				and ( mVectorOffset == other.mVectorOffset );
		}

		/**
		 * Compare iterators for inequality.
		 * @param other Const reference to the iterator to compare against for inequality.
		 * @return True is returned if this and other are not equal.
		 */
		bool operator!=(
			const iterator& other )
		{
			return not this->operator==( other );
		}

		/**
		 * Pointer access to the FastaSequence.
		 * @return A pointer to the underlying FastaSequence object.
		 */
		pointer operator->()
		{
			if ( mIsValid )
			{
				return &mMapIterator.operator*().second[ mVectorOffset ];
			}

			return static_cast< pointer >( nullptr );
		}

		/**
		 * Reference acess to the FastaSequence.
		 * @return A reference to the underlying FastaSequence object.
		 */
		reference operator*()
		{
			if ( mIsValid )
			{
				return mMapIterator.operator*().second[ mVectorOffset ];
			}

			return *static_cast< pointer >( nullptr );
		}

		/**
		 * Post-increment operator.
		 * @return Return the iterator prior to incrementing.
		 */
		iterator operator++( int )
		{
			iterator prior( *this );
			this->operator++();
			return prior;
		}

		/**
		 * Pre-increment operator.
		 * @return Reference to this iterator instance is returned.
		 */
		iterator& operator++()
		{
			if ( mIsValid )
			{
				mMapIterator.operator++();
			}

			return *this;
		}

		/**
		 * Swap the contents of this iterator with other.
		 * @param other Reference to the iterator instance to swap values with.
		 */
		void swap(
			iterator& other )
		{
			std::swap( mIsValid, other.mIsValid );
			std::swap( mMapIterator, other.mMapIterator );
			std::swap( mVectorOffset, other.mVectorOffset );
		}
	};

	/**
	 * Class for iterating over the container in an immutable manner.
	 */
	class const_iterator
	{
	private:
		friend class FastaFile;

		bool mIsValid;
		IdentifierSequenceMapType::const_iterator mConstMapIterator;
		size_t mVectorOffset;

		const_iterator(
			IdentifierSequenceMapType::const_iterator constMapIterator )
		{
			mIsValid = true;
			mConstMapIterator = constMapIterator;
			mVectorOffset = 0;
		}

		void _copyAssign(
			const const_iterator& other )
		{
			mIsValid = other.mIsValid;
			mConstMapIterator = other.mConstMapIterator;
			mVectorOffset = other.mVectorOffset;
		}

		void _moveAssign(
			const_iterator&& other )
		{
			mIsValid = std::exchange( other.mIsValid, false );
			mConstMapIterator = std::move( other.mConstMapIterator );
			mVectorOffset = std::exchange( other.mVectorOffset, 0 );
		}

	public:
		using iterator_category = std::forward_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = FastaSequence;
		using pointer           = const FastaSequence*;
		using reference         = const FastaSequence&;

		/**
		 * Default constructor to an empty const_iterator.
		 */
		const_iterator()
		{
			mIsValid = false;
			mVectorOffset = 0;
		}

		/**
		 * Copy constructor.
		 * @param other Const reference to const_iterator to copy to this instance.
		 */
		const_iterator(
			const const_iterator& other )
		{
			_copyAssign( other );
		}

		/**
		 * Move constructor.
		 * @param other R-Value to const_iterator to move to this instance.
		 */
		const_iterator(
			const_iterator&& other )
		{
			_moveAssign( std::move( other ) );
		}

		/**
		 * Copy assignment operator.
		 * @param other Const reference to const_iterator to copy to this instance.
		 * @return Reference to this const_iterator instance is returned.
		 */
		const_iterator& operator=(
			const const_iterator& other )
		{
			if ( this != &other )
			{
				_copyAssign( other );
			}

			return *this;
		}

		/**
		 * Move assignment operator.
		 * @param other R-Value to const_iterator to move to this instance.
		 * @return Reference to this const_iterator instance is returned.
		 */
		const_iterator& operator=(
			const_iterator&& other )
		{
			if ( this != &other )
			{
				_moveAssign( std::move( other ) );
			}

			return *this;
		}

		/**
		 * Compare const_iterators for equality.
		 * @param other Const reference to the const_iterator to compare against for equality.
		 * @return True is returned if this and other are equal.
		 */
		bool operator==(
			const const_iterator& other )
		{
			return mIsValid
				and ( mIsValid == other.mIsValid )
				and ( mConstMapIterator == other.mConstMapIterator )
				and ( mVectorOffset == other.mVectorOffset );
		}

		/**
		 * Compare const_iterators for inequality.
		 * @param other Const reference to the const_iterator to compare against for inequality.
		 * @return True is returned if this and other are not equal.
		 */
		bool operator!=(
			const const_iterator& other )
		{
			return not this->operator==( other );
		}

		/**
		 * Const pointer access to the FastaSequence.
		 * @return A const pointer to the underlying FastaSequence object.
		 */
		pointer operator->()
		{
			if ( mIsValid )
			{
				return &mConstMapIterator.operator*().second[ mVectorOffset ];
			}

			return static_cast< pointer >( nullptr );
		}

		/**
		 * Const reference acess to the FastaSequence.
		 * @return A const reference to the underlying FastaSequence object.
		 */
		reference operator*()
		{
			if ( mIsValid )
			{
				return mConstMapIterator.operator*().second[ mVectorOffset ];
			}

			return *static_cast< pointer >( nullptr );
		}

		/**
		 * Post-increment operator.
		 * @return Return the const_iterator prior to incrementing.
		 */
		const_iterator operator++( int )
		{
			const_iterator prior( *this );
			this->operator++();
			return prior;
		}

		/**
		 * Pre-increment operator.
		 * @return Reference to this const_iterator instance is returned.
		 */
		const_iterator& operator++()
		{
			if ( mIsValid )
			{
				mConstMapIterator.operator++();
			}

			return *this;
		}

		/**
		 * Swap the contents of this const_iterator with other.
		 * @param other Reference to the const_iterator instance to swap values with.
		 */
		void swap(
			const_iterator& other )
		{
			std::swap( mIsValid, other.mIsValid );
			std::swap( mConstMapIterator, other.mConstMapIterator );
			std::swap( mVectorOffset, other.mVectorOffset );
		}
	};

	/**
	 * Default constructor.
	 * @param filename The name of the file to load into this instance. [default: ""]
	 * @param allowDuplicates Flag to allow or disallow duplicate identifiers. [default: true]
	 */
	FastaFile(
		const std::string& filename = std::string(),
		bool allowDuplicates = true )
	{
		mIsBareSequence = false;

		if ( filename.empty() )
		{
			mDuplicateIdentifiersAllowed = allowDuplicates;
		}
		else
		{
			this->readFile( filename, allowDuplicates );
		}
	}

	/**
	 * Copy constructor.
	 * @param other Const reference to the FastaFile to copy to this instance.
	 */
	FastaFile(
		const FastaFile& other )
	{
		_copyAssign( other );
	}

	/**
	 * Move constructor.
	 * @param other R-Value to the FastaFile to move to this instance.
	 */
	FastaFile(
		FastaFile&& other )
	{
		_moveAssign( std::move( other ) );
	}

	/**
	 * Add the sequence to the container.
	 * @param Const reference to the sequence to add to the container.
	 * @return If duplicates are allowed, then this will always return 1.
	 *         If duplicates are not allowed, then 1 will only be returned if the
	 *         identifier for the sequence is not already taken.
	 */
	size_t addSequence(
		const FastaSequence& sequence )
	{
		return this->addSequences( std::vector< FastaSequence >{ sequence } );
	}

	/**
	 * Add the sequence to the container.
	 * @param Const reference to the sequence to add to the container.
	 * @return If duplicates are allowed, then this will always return 1.
	 *         If duplicates are not allowed, then 1 will only be returned if the
	 *         identifier for the sequence is not already taken.
	 */
	size_t addSequence(
		FastaSequence&& sequence )
	{
		return this->addSequences( std::vector< FastaSequence >{ std::move( sequence ) } );
	}

	/**
	 * Add a vector of sequences to the container.
	 * @param Const reference to the vector of sequences to add to the container.
	 * @return The number of sequences from the vector added to the container from the
	 *         vector is returned. If duplicates are not allowed, then the number returned
	 *         may be less than the size of the vector.
	 */
	size_t addSequences(
		const std::vector< FastaSequence >& sequences )
	{
		size_t numberAdded( 0 );

		if ( mDuplicateIdentifiersAllowed )
		{
			for ( const auto& sequence : sequences )
			{
				mIdentifierSequenceMap[ sequence.identifier() ].push_back( sequence );
			}

			numberAdded = sequences.size();
		}
		else
		{
			for ( const auto& sequence : sequences )
			{
				if ( mIdentifierSequenceMap.end() == mIdentifierSequenceMap.find( sequence.identifier() ) )
				{
					mIdentifierSequenceMap[ sequence.identifier() ].push_back( sequence );
					++numberAdded;
				}
			}
		}

		return numberAdded;
	}

	/**
	 * Flag that the container is or is not allowed duplicate identifiers.
	 * By default, the container allows duplicate identifiers. Should false be
	 * passed into this method, then a scan is performed to removed any
	 * sequences that map to a shared identifier excluding the first added to this container.
	 * @param allow Flag whether or not duplicate identifiers are allowed.
	 */
	void allowDuplicateIdentifiers(
		bool allow = true )
	{
		mDuplicateIdentifiersAllowed = allow;

		if ( not mDuplicateIdentifiersAllowed )
		{
			for ( auto& sequenceVector : mIdentifierSequenceMap )
			{
				sequenceVector.second.resize( 1 );
			}
		}
	}

	/**
	 * Get the vector of FastaSequences associated with the given identifier.
	 * @param identifier The identifier to the associated sequences.
	 * @return Const reference to the vector of FastaSequences associated with the given identifier.
	 * @throw std::out_of_range is thrown if no such identifier is present in the container.
	 */
	const std::vector< FastaSequence >& at(
		const std::string& identifier ) const
	{
		return mIdentifierSequenceMap.at( identifier );
	}

	/**
	 * Get a const_iterator to the beginning of the container.
	 * @return A const_iterator to the beginning of the container is returned.
	 */
	const_iterator begin() const
	{
		return const_iterator( mIdentifierSequenceMap.begin() );
	}

	/**
	 * Get an iterator to the beginning of the container.
	 * @return An iterator to the beginning of the container is returned.
	 */
	iterator begin()
	{
		return iterator( mIdentifierSequenceMap.begin() );
	}

	/**
	 * Get a const_iterator to the beginning of the container.
	 * @return A const_iterator to the beginning of the container is returned.
	 */
	const_iterator cbegin() const
	{
		return const_iterator( mIdentifierSequenceMap.begin() );
	}

	/**
	 * Get a const_iterator to the end of the container.
	 * @return A const_iterator to the end of the container is returned.
	 */
	const_iterator cend() const
	{
		return const_iterator( mIdentifierSequenceMap.end() );
	}

	/**
	 * Get a const_iterator to the end of the container.
	 * @return A const_iterator to the end of the container is returned.
	 */
	const_iterator end() const
	{
		return const_iterator( mIdentifierSequenceMap.end() );
	}

	/**
	 * Get an iterator to the end of the container.
	 * @return An iterator to the end of the container is returned.
	 */
	iterator end()
	{
		return iterator( mIdentifierSequenceMap.end() );
	}

	/**
	 * Retrieve the list of identifiers present within the container.
	 * @return A const reference to the vector of identifiers present.
	 */
	const std::set< std::string >& getIdentifiers() const
	{
		return mIdentifiersSet;
	}

	/**
	 * Check for the presence of an identifier in the container.
	 * @param identifier Const reference to the identifier to check for.
	 * @return True is returned if the identifier is present.
	 */
	bool hasIdentifier(
		const std::string& identifier ) const
	{
		return mIdentifiersSet.end() != mIdentifiersSet.find( identifier );
	}

	/**
	 * Returns true if the file is a bare sequence.
	 * @return True is returned if the file is a bare sequence.
	 */
	bool isBareSequence() const
	{
		return mIsBareSequence;
	}

	/**
	 * Copy assignment operator.
	 * @param other Const reference to the FastaFile to copy to this instance.
	 * @return Reference to this FastaFile instance is returned.
	 */
	FastaFile& operator=(
		const FastaFile& other )
	{
		if ( this != &other )
		{
			_copyAssign( other );
		}

		return *this;
	}

	/**
	 * Move assignment operator.
	 * @param other R-Value to the FastaFile to move to this instance.
	 * @return Reference to this FastaFile instance is returned.
	 */
	FastaFile& operator=(
		FastaFile&& other )
	{
		if ( this != &other )
		{
			_moveAssign( std::move( other ) );
		}

		return *this;
	}

	/**
	 * Read in a FastA file into this FastaFile instance. If {@param allowDuplicates} is set
	 * to false and there are duplicates present in the file, then only the first sequence is selected.
	 * @param filename The name of the file to load into this instance.
	 * @param allowDuplicates Flag to allow or disallow duplicate identifiers in the source file. [default: true]
	 * @return Zero is returned upon success, else an error code is returned.
	 */
	int readFile(
		const std::string& filename,
		bool allowDuplicates = true )
	{
		std::string line, sequenceString;
		std::ifstream inputFile;
		FastaSequence sequence;

		inputFile.open( filename, std::ios::in );

		while ( std::getline( inputFile, line ) )
		{
			if ( '>' == line[ 0 ] )
			{
				if ( 0 < sequence.identifier().length() )
				{
					sequence.setSequence( sequenceString );
					mIdentifiersSet.insert( sequence.identifier() );
					mIdentifierSequenceMap[ sequence.identifier() ].push_back( sequence );
				}

				sequence.setIdentifier( line );
				sequenceString.clear();
			}
			else
			{
				sequenceString.append( line );
			}
		}

		sequence.setSequence( sequenceString );
		mIdentifiersSet.insert( sequence.identifier() );
		mIdentifierSequenceMap[ sequence.identifier() ].push_back( sequence );

		inputFile.close();

		this->allowDuplicateIdentifiers( allowDuplicates );

		return 0;
	}

	/**
	 * Write the contents of this container out to file.
	 * @param filename Name of the file to write to.
	 * @param lineLength Length of each sequence line. If {@param lineLength} is zero, then
	 *                   the entire sequence will be on a single line. [default: 80]
	 * @return Return zero upon success, else an error code is returned.
	 */
	int writeFile(
		const std::string& filename,
		size_t lineLength = 80 ) const
	{
		std::ofstream outputFile;
		outputFile.open( filename, std::ios::out | std::ios::trunc );

		if ( 0 == lineLength )
		{
			lineLength = static_cast< size_t >( -1 );
		}

		for ( const auto& sequenceVector : mIdentifierSequenceMap )
		{
			for ( const auto& sequence : sequenceVector.second )
			{
				outputFile << ">" << sequence.identifier() << std::endl;

				for ( size_t offset( 0 ); offset < sequence.length(); offset += lineLength )
				{
					outputFile << sequence.sequence().substr( offset, lineLength ) << std::endl;
				}
			}
		}

		outputFile.close();

		return 0;
	}
};
