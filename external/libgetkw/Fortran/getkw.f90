! getkw -- a simple input parser
!
! Copyright (C) 2006  Jonas Juselius
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
! Wrapper module for the parser. Uses the wonderful f90 overloading caps.
! 

module getkw_class
    use kinds_m
	use teletype_m
    implicit none
	
	public new_getkw, del_getkw, getkw, getkw_ptr, setkw, save_keys
	public push_section, pop_section, has_keyword
	public keyword_is_set, section_is_set
	public addkw, add_section, delkw, del_section
	public getsect, set_verbose, set_strict
	public kword_t, section_t, getkw_t
	public memtrack, typtrack
	public LINELEN
	public print_tree
	private

	interface getkw
		module procedure getkw_ival
		module procedure getkw_dval
		module procedure getkw_lval
		module procedure getkw_ivec
		module procedure getkw_dvec
		module procedure getkw_lvec
		module procedure getkw_string
		module procedure getkw_str
		module procedure getkw_ival_ptr
		module procedure getkw_ivec_ptr
		module procedure getkw_dval_ptr
		module procedure getkw_dvec_ptr
		module procedure getkw_lval_ptr
		module procedure getkw_lvec_ptr
		module procedure getkw_str_ptr
	end interface

	interface getkw_ptr
		module procedure getkw_str_ref
		module procedure getkw_lvec_ref
		module procedure getkw_dvec_ref
		module procedure getkw_ivec_ref
	end interface

	interface setkw
		module procedure setkw_ival
		module procedure setkw_ivec
		module procedure setkw_dval
		module procedure setkw_dvec
		module procedure setkw_lval
		module procedure setkw_lvec
		module procedure setkw_str
	end interface

	interface addkw
		module procedure addkw_ival
		module procedure addkw_ivec
		module procedure addkw_dval
		module procedure addkw_dvec
		module procedure addkw_lval
		module procedure addkw_lvec
		module procedure addkw_string
		module procedure addkw_strvec
	end interface
	
    integer, parameter :: MAXID=32
    integer, parameter :: LINELEN=132
    integer, parameter :: MAX_STACK=10

    integer, parameter :: KW_NONE=-1
    integer, parameter :: KW_INT=1
    integer, parameter :: KW_DBL=2
    integer, parameter :: KW_BOOL=3
    integer, parameter :: KW_STR=4
    integer, parameter :: KW_IVEC=5
    integer, parameter :: KW_DVEC=6
    integer, parameter :: KW_LVEC=7

	! this really ought to be a union of sorts...
    type kword_t
        private
		real(DP) :: dval
		integer(SP) :: ival
		logical :: bool, set=.false.
		character(LINELEN), dimension(:), pointer :: str
		real(DP), dimension(:), pointer :: dvec
		integer(SP), dimension(:), pointer :: ivec
		logical, dimension(:), pointer :: lvec
		type(kword_t), pointer :: next
    end type

	type keyword_t
		private
        logical :: set
		integer(SP) :: typ
        character(MAXID) :: id
		integer :: nkeys
		type(kword_t), pointer :: key
		type(keyword_t), pointer :: next
	end type

    type section_t
        private
        logical :: set
		integer(SP) :: nsect, nkw
        character(MAXID) :: id
        type(keyword_t), pointer :: arg
        type(keyword_t), pointer :: kw
		type(section_t), pointer :: sects
		type(section_t), pointer :: next
    end type

	type sectstack_t
		type(section_t), pointer :: sect
	end type

    type getkw_t
        private
        type(section_t), pointer :: main
		type(section_t), pointer :: active
		type(sectstack_t), dimension(MAX_STACK)  :: stack
		integer :: stidx
		integer :: lun=5
!        type(memstat_t) :: mem
    end type

	type line_t
		integer :: len
		character, dimension(:), pointer :: foo
	end type

	type(getkw_t), pointer :: main
	type(section_t), pointer :: active
	integer :: memtrack=0, typtrack=0
	integer :: lun=5

	logical :: verbose=.true.
	logical :: strict=.false.

contains

	subroutine globals(self)
		type(getkw_t), target :: self

		lun=self%lun
		main=>self
		active=>self%active
	end subroutine

	subroutine new_keyword(key)
		type(keyword_t), pointer :: key

		allocate(key)
		typtrack=typtrack+1
		nullify(key%next)
	end subroutine

	subroutine new_kword(key)
		type(kword_t), pointer :: key

		allocate(key)
		typtrack=typtrack+1
		nullify(key%next)
		nullify(key%ivec)
		nullify(key%dvec)
		nullify(key%lvec)
		nullify(key%str)
	end subroutine

	subroutine new_section(sect)
		type(section_t), pointer :: sect

		allocate(sect)
		typtrack=typtrack+1
		nullify(sect%arg)
		nullify(sect%kw)
		nullify(sect%sects)
		nullify(sect%next)
	end subroutine

	subroutine delkw(self, path)
		type(getkw_t) :: self
		character(*), intent(in) :: path

		logical :: ok
		integer(SP) :: idx
		type(section_t), pointer :: sect
		type(keyword_t), pointer :: kw, kw2

		idx=index(path, '.', back=.true.)
		if (idx > 0) then
			ok=find_sect(self%active, path(:idx-1), sect)
			if (.not.ok) then
				call msg_error('delkw: no such key: ' // trim(path))
				return
			end if
		else
			sect=>self%active
		end if

		ok=findkw(self%active, path, kw, kw2)
		if (.not.ok) then
			call msg_error('delkw: no such key: ' // trim(path))
			return
		end if
		kw=>del_keyword(kw)
		if (associated(kw2)) then
			kw2%next=>kw
		else
			sect%kw=>kw
		end if
		sect%nkw=sect%nkw-1
	end subroutine

	subroutine del_section(self, path)
		type(getkw_t) :: self
		character(*), intent(in) :: path

		logical :: ok
		integer(SP) :: idx
		type(section_t), pointer :: ptr, sect, sect2

		idx=index(path, '.', back=.true.)
		if (idx > 0) then
			ok=find_sect(self%active, path(:idx-1), ptr)
			if (.not.ok) then
				call msg_error('delkw: no such section: ' // trim(path))
				return
			end if
		else
			ptr=>self%active
		end if

		ok=find_sect(ptr, path, sect, sect2)
		if (.not.ok) then
			call msg_error('delkw: no such key: ' // trim(path))
			return
		end if
		sect=>del_section_t(sect)
		if (associated(sect2)) then
			sect2%next=>sect
		else
			ptr%sects=>sect
		end if
		sect%nsect=sect%nsect-1
	end subroutine

	function del_keyword(key) result(next)
		type(keyword_t), pointer :: key, next
		
		next=>key%next
		if (associated(key%key))  then
			call del_kword(key%key)
		else
			call msg_warn('Cannot delete keyword: not associated!')
		end if

		deallocate(key)
		typtrack=typtrack-1
	end function

	recursive subroutine del_kword(key) 
		type(kword_t), pointer :: key
		
		if (associated(key%str)) then
			deallocate(key%str)
			memtrack=memtrack-1
		end if
		if (associated(key%ivec)) then
			deallocate(key%ivec)
			memtrack=memtrack-1
		end if
		if (associated(key%dvec)) then
			deallocate(key%dvec)
			memtrack=memtrack-1
		end if
		if (associated(key%lvec)) then
			deallocate(key%lvec)
			memtrack=memtrack-1
		end if
		
		do while (associated(key%next))
			call del_kword(key%next)
		end do

		deallocate(key)
		typtrack=typtrack-1
	end subroutine

	recursive function del_section_t(sect) result(next)
		type(section_t), pointer :: sect
		type(section_t), pointer :: next
		
		type(keyword_t), pointer :: kw 

		if (associated(sect%arg)) then
			kw=>del_keyword(sect%arg)
		end if

		if (associated(sect%kw)) then
			kw=>sect%kw
			do while (associated(kw))
				kw=>del_keyword(kw)
			end do
		end if

		if (associated(sect%sects)) then
			next=>sect%sects
			do while (associated(next))
				next=>del_section_t(next)
			end do
		end if

		next=>sect%next
		deallocate(sect)
		typtrack=typtrack-1
	end function

    subroutine new_getkw(self, ilun)
		type(getkw_t), target :: self
		integer(SP), optional :: ilun
		
		if (present(ilun)) then
			self%lun=ilun
		end if

		call new_section(self%main)
		call read_section(self%main)
        self%active=>self%main

!        call print_section(self%main)
!        call globals(self)
    end subroutine
	

	subroutine del_getkw(self)
		type(getkw_t) :: self
		type(section_t), pointer :: main

		if (.not.associated(self%main)) then 
			call msg_error('del_getkw: main section not associated!')
			stop
		end if
		main=>self%main
		do while (associated(main))
			main=>del_section_t(main)
		end do
		nullify(main)
	end subroutine

!::::::::::::::::::::  getkw  :::::::::::::::::::::::::::

	subroutine getkw_err(str)
		character(*), intent(in) :: str

		if (strict) then
			call msg_error('getkw: no such key: ' // trim(str))
			stop
		else if (verbose) then
			call msg_warn('getkw: no such key: ' // trim(str))
		end if
	end subroutine

	subroutine getkw_errundef(str)
		character(*), intent(in) :: str

		if (strict) then
			call msg_error('getkw: key undefined: ' // trim(str))
			stop
		else if (verbose) then
			call msg_warn('getkw: key undefined: ' // trim(str))
		end if
	end subroutine

	subroutine getkw_typerr(kw)
		type(keyword_t) :: kw

		if (strict) then
			call msg_error('getkw: invalid kw type: '//trim(kw%id)// ' -> ' // &
			xstr(kw%typ))
			stop
		else if (verbose) then
			call msg_warn('getkw: invalid kw type: '//trim(kw%id)// ' -> ' // &
			xstr(kw%typ))
		end if
	end subroutine

	subroutine getkw_ival(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr

		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_INT) then
				call getkw_typerr(ptr)
			else
				val=ptr%key%ival
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

	subroutine getkw_ivec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP), dimension(:) :: val
		
		integer(SP), dimension(:), pointer :: ptr

		nullify(ptr)
		call getkw_ivec_ref(self, path, ptr)
		if (associated(ptr)) val=ptr
	end subroutine

	subroutine getkw_ivec_ref(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP), dimension(:), pointer :: val
		
		logical :: ok
		type(keyword_t), pointer :: ptr

		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_IVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%ivec
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

	subroutine getkw_dvec_ref(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_DVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%dvec
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

	subroutine getkw_dvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), dimension(:) :: val
		
		real(DP), dimension(:), pointer :: ptr
		nullify(ptr)

		call getkw_dvec_ref(self, path, ptr)
		if (associated(ptr)) val=ptr
	end subroutine

	subroutine getkw_lvec_ref(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_LVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%lvec
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

	subroutine getkw_lvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, dimension(:) :: val
		
		logical, dimension(:), pointer :: ptr
		nullify(ptr)

		call getkw_lvec_ref(self, path, ptr)
		if (associated(ptr)) val=ptr
	end subroutine

	subroutine getkw_str_ref(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		character(*), dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_STR) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%str
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

	subroutine getkw_str(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		character(LINELEN), dimension(:), pointer :: val

		call getkw_str_ref(self, path, val)
!        if (associated(ptr)) val=ptr
	end subroutine

	subroutine getkw_string(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		character(*), intent(out) :: val

		character(LINELEN), dimension(:), pointer :: ptr
		nullify(ptr)
		call getkw_str_ref(self, path, ptr)
		if (associated(ptr)) val=ptr(1)
	end subroutine
	
	subroutine getkw_dval(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), intent(out) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_DBL) then
				call getkw_typerr(ptr)
			else
				val=ptr%key%dval
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

	subroutine getkw_lval(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, intent(out) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_BOOL) then
				call getkw_typerr(ptr)
			else
				val=ptr%key%bool
			end if
		else
			call getkw_err(path)
		end if
	end subroutine

!::::::::::::::::::::  getkw_ptr  :::::::::::::::::::::::::::

	subroutine getkw_ival_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		integer(SP), intent(out) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_INT) then
				call getkw_typerr(ptr)
			else
				val=ptr%key%ival
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine


	subroutine getkw_ivec_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		integer(SP), dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_IVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%ivec
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine

	subroutine getkw_dvec_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		real(DP), dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_DVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%dvec
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine

	subroutine getkw_lvec_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		logical, dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_LVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%lvec
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine

	subroutine getkw_str_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		character(*), dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_LVEC) then
				call getkw_typerr(ptr)
			else
				val=>ptr%key%str
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine
	
	subroutine getkw_dval_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		real(DP), intent(out) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_DBL) then
				call getkw_typerr(ptr)
			else
				val=ptr%key%dval
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine

	subroutine getkw_lval_ptr(self, sect, val)
		type(getkw_t), target :: self
		type(section_t), pointer :: sect
		logical, intent(out) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(sect, ' ', ptr)
		if (ok) then
            if (.not.ptr%key%set) then
!                call getkw_errundef(ptr%id)
				return
			end if
			if (ptr%typ	 /= KW_BOOL) then
				call getkw_typerr(ptr)
			else
				val=ptr%key%bool
			end if
		else
			call getkw_err(sect%id)
		end if
	end subroutine


!::::::::::::::::::::  addkw  ::::::::::::::::::::::::::::
	subroutine addkw_err(str)
		character(*), intent(in) :: str
		if (strict) then
			call msg_error('addkw: unknown section: ' // trim(str)) 
			stop
		else if (verbose) then
			call msg_warn('addkw: unknown section: ' // trim(str)) 
		end if
	end subroutine

	function insert_keyword(sect, path) result(kw)
		type(section_t), pointer :: sect
		character(*), intent(in) :: path
		type(keyword_t), pointer :: kw

		logical :: ok, arg
		integer(SP) :: idx
		type(section_t), pointer :: ptr
		nullify(ptr)

		ok=.false.
		arg=.false.
		idx=index(path, '.', back=.true.)
		if (idx > 0) then
			ok=find_sect(sect, path(:idx-1), ptr)
		else
			ptr=>sect
			ok=.true.
		end if

		if (.not.ok) return

		ok=findkw(ptr, path(idx+1:), kw)
		if (ok) then
			if (.not.associated(kw)) then ! section arg
				ok=find_sect(sect, path, ptr)
				arg=.true.
			else
				call msg_error('insert_keyword: key already defined: '// &
				trim(path))
				ok=.false.
				stop  ! is this right???
			end if
		end if

		if (arg) then
			call new_keyword(ptr%arg)
			kw=>ptr%arg
		else
			if (.not.associated(ptr%kw)) then
				call new_keyword(ptr%kw)
				kw=>ptr%kw
			else
				kw=>last_keyword(ptr%kw)
				call new_keyword(kw%next)
				kw=>kw%next
			end if
		end if

		if (idx > 0) then
			kw%id=path(idx+1:)
		else
			kw%id=path
		end if
		kw%nkeys=1
		kw%set=.true.
		kw%typ=-1
		call new_kword(kw%key)
		if (.not.arg) then
			ptr%nkw=ptr%nkw+1
		end if
	end function

	subroutine addkw_ival(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
			kw%key%ival=val
			kw%typ=KW_INT
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_ivec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP), dimension(:), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%typ=KW_IVEC
				allocate(kw%key%ivec(size(val)))
				memtrack=memtrack+1
				kw%key%ivec=val
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_dval(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%key%dval=val
				kw%typ=KW_DBL
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_dvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), dimension(:), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%typ=KW_DVEC
				allocate(kw%key%dvec(size(val)))
				memtrack=memtrack+1
				kw%key%dvec=val
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_lval(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%key%bool=val
				kw%typ=KW_BOOL
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_lvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, dimension(:), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%typ=KW_LVEC
				allocate(kw%key%lvec(size(val)))
				memtrack=memtrack+1
				kw%key%lvec=val
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_string(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		character(*), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%typ=KW_STR
				allocate(kw%key%str(1))
				memtrack=memtrack+1
				kw%key%str=val
		else
			call addkw_err(path)
		end if
	end subroutine

	subroutine addkw_strvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		character(LINELEN), dimension(:), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: kw

        kw=>insert_keyword(self%active, path)
		
		if (associated(kw)) then
				kw%typ=KW_STR
				allocate(kw%key%str(size(val)))
				memtrack=memtrack+1
				kw%key%str=val
		else
			call addkw_err(path)
		end if
	end subroutine

!::::::::::::::::::::  add_section  ::::::::::::::::::::::::::::
	function insert_section(sect, path) result(new)
		type(section_t), pointer :: sect
		character(*), intent(in) :: path
		type(section_t), pointer :: new

		logical :: ok
		integer(SP) :: idx
		type(section_t), pointer :: ptr
		nullify(ptr)

		ok=.false.
		idx=index(path, '.', back=.true.)
		if (idx > 0) then
			ok=find_sect(sect, path(:idx-1), ptr)
		else
			ptr=>sect
			ok=.true.
		end if

		if (ok) then
			if (.not.associated(ptr%sects)) then
				call new_section(ptr%sects)
				new=>ptr%sects
			else
				new=>last_section(ptr%sects)
				call new_section(new%next)
				new=>new%next
			end if
			if (idx > 0) then
				new%id=path(idx+1:)
			else
				new%id=path
			end if
			new%set=.true.
			new%nkw=0
			new%nsect=0
			ptr%nsect=ptr%nsect+1
		else
			nullify(new)
		end if
	end function

	subroutine add_section(self, path)
		type(getkw_t) :: self
		character(*), intent(in) :: path

		type(section_t), pointer :: new

		new=>insert_section(self%active, path)
		if (.not.associated(new)) then
			call addkw_err(path)
		end if
	end subroutine

!::::::::::::::::::::  setkw  ::::::::::::::::::::::::::::
	subroutine setkw_err(str)
		character(*), intent(in) :: str
		if (strict) then
			call msg_error('setkw: unknown key: ' // trim(str)) 
			stop
		else if (verbose) then
			call msg_warn('setkw: unknown key: ' // trim(str)) 
		end if
	end subroutine

	subroutine setkw_ival(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_INT) then
				call getkw_typerr(ptr)
			else
				ptr%key%ival=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

	subroutine setkw_ivec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		integer(SP), dimension(:) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_IVEC) then
				call getkw_typerr(ptr)
			else
				if (associated(ptr%key%ivec)) then
					deallocate(ptr%key%ivec)
				end if
				allocate(ptr%key%ivec(size(val)))
				ptr%key%ivec=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

	subroutine setkw_dvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), dimension(:) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_DVEC) then
				call getkw_typerr(ptr)
			else
				if (associated(ptr%key%dvec)) then
					deallocate(ptr%key%dvec)
				end if
				allocate(ptr%key%dvec(size(val)))
				ptr%key%dvec=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

	subroutine setkw_dval(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		real(DP), intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_DBL) then
				call getkw_typerr(ptr)
			else
				ptr%key%dval=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

	subroutine setkw_lvec(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_LVEC) then
				call getkw_typerr(ptr)
			else
				if (associated(ptr%key%lvec)) then
					deallocate(ptr%key%lvec)
				end if
				allocate(ptr%key%lvec(size(val)))
				ptr%key%lvec=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

	subroutine setkw_lval(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		logical, intent(in) :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_BOOL) then
				call getkw_typerr(ptr)
			else
				ptr%key%bool=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

	subroutine setkw_str(self, path, val)
		type(getkw_t), target :: self
		character(*), intent(in) :: path
		character(*), dimension(:), pointer :: val
		
		logical :: ok
		
		type(keyword_t), pointer :: ptr
		nullify(ptr)
		call globals(self)

		ok=findkw(self%active, path, ptr)
		if (ok) then
			if (ptr%typ	 /= KW_LVEC) then
				call getkw_typerr(ptr)
			else
				if (associated(ptr%key%str)) then
					deallocate(ptr%key%str)
				end if
				allocate(ptr%key%str(size(val)))
				ptr%key%str=val
				ptr%set=.true.
			end if
		else
			call setkw_err(path)
		end if
	end subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	subroutine save_keys(self, lun)
		type(getkw_t) :: self
		integer(SP), intent(in) :: lun

		type(section_t), pointer :: foo	

		call globals(self)

		! for all sections: print section
		! for all keys in section: print key
		call write_section(self%active, lun)
	end subroutine

	subroutine write_key(kw, lun) 
		type(keyword_t), pointer :: kw
		integer(SP) :: lun

		integer(SP) :: i, nstr

		select case (kw%typ)
			case(KW_INT)
				write(lun,*) 'INT ', trim(kw%id), 1, kw%set
				write(lun,*) kw%key%ival
			case(KW_DBL)
				write(lun,*) 'DBL ', trim(kw%id), 1, kw%set
				write(lun,*) kw%key%dval
			case(KW_BOOL)
				write(lun,*) 'BOOL ', trim(kw%id), 1, kw%set
				write(lun,*) kw%key%bool
			case(KW_IVEC)
				write(lun,*) 'INT ', trim(kw%id), size(kw%key%ivec), kw%set
				write(lun,*) kw%key%ivec
			case(KW_DVEC)
				write(lun,*) 'DBL ', trim(kw%id), size(kw%key%dvec), kw%set
				write(lun,*) kw%key%dvec
			case(KW_LVEC)
				write(lun,*) 'BOOL ', trim(kw%id), size(kw%key%lvec), kw%set
				write(lun,*) kw%key%lvec
			case(KW_STR)
				if (associated(kw%key%str)) then
					nstr=size(kw%key%str)
					write(lun,*) 'STR ', trim(kw%id), nstr, kw%set
					do i=1,nstr
						write(lun,*) trim(kw%key%str(i))
					end do
				else
					write(lun,*) 'STR ', trim(kw%id), 0, kw%set
				end if
			case(KW_NONE)
				call msg_warn('key undefined')
			case default
				call msg_warn('dang! no such key type' // xstr(kw%typ))
		end select
	end subroutine

	recursive subroutine write_section(sect, lun)
		type(section_t) :: sect
		integer(SP) :: lun

		logical :: arg
		integer(SP) :: nkw, i
		type(keyword_t), pointer :: kw
		type(section_t), pointer :: ptr
		nullify(ptr)

		write(lun,*) 'SECT ', trim(sect%id), sect%nsect, sect%set
		write(lun,*) 'ARG ', associated(sect%arg), ' KW ', sect%nkw

		if (associated(sect%arg)) then 
			call write_key(sect%arg, lun) 
		end if
		
		kw=>sect%kw
		do while (associated(kw))
			call write_key(kw, lun)
			kw=>next_keyword(kw)
		end do

		ptr=>sect%sects
		do while (associated(ptr))
			call write_section(ptr, lun)
			ptr=>next_section(ptr)
		end do
	end subroutine

	recursive function findkw(cur, path, ptr, ptr2) result(ok)
		type(section_t) :: cur
		character(*), intent(in) :: path
		type(keyword_t), pointer :: ptr
		type(keyword_t), optional, pointer :: ptr2
		logical :: ok

		type(section_t), pointer :: next
		integer(SP) :: index, idx, i, nkw
		
		nullify(ptr)

		ok=.false.
			
		idx=index(path, '.')
		if (idx > 0) then
			next=>cur%sects
			do while (associated(next))
				if (path(:idx-1) == next%id) then
					ok=findkw(next, path(idx+1:), ptr)
				end if
				if (ok) exit
				next=>next_section(next)
			end do
		else
			ptr=>cur%kw
			if (present(ptr2)) nullify(ptr2)
			do while (associated(ptr))
				if (path == ptr%id) then
					ok=.true.
					return
				end if
				if (present(ptr2)) ptr2=>ptr
				ptr=>ptr%next
			end do
			! hmm, let's see if it's a section arg then...
			ptr=>cur%arg
			next=>cur%sects
			do while (associated(next))
				if (path == next%id) then
					ptr=>next%arg
					if (associated(next%arg)) then
						ptr=>next%arg
						ok=.true.
					else
						call msg_warn('no argument for section: '//path)
						ok=.false.
					end if
					return
				else
					next=>next_section(next)
				end if
			end do
			nullify(ptr)
		end if
	end function

	subroutine getsect(self, path, ptr)
		type(getkw_t) :: self
		character(*), intent(in) :: path
		type(section_t), pointer :: ptr

		logical :: ok

		call globals(self)
		ok=find_sect(self%active, path, ptr)
		if (.not.ok) then
			nullify(ptr)
		end if
	end subroutine

	recursive function find_sect(cur, path, ptr, ptr2) result(ok)
		type(section_t) :: cur
		character(*), intent(in) :: path
		type(section_t), pointer :: ptr, new
		type(section_t), optional, pointer :: ptr2
		logical :: ok

		integer(SP) :: index, idx, i

		ok=.false.
		idx=index(path, '.')
		ptr=>cur%sects
		if (idx > 0) then
			do while (associated(ptr))
				if (path(:idx-1) == ptr%id) then
					ok=find_sect(ptr, path(idx+1:), new)
				end if
				if (ok) then
					ptr=>new
					exit
				end if
				ptr=>next_section(ptr)
			end do
		else
			if (present(ptr2)) nullify(ptr2)
			do while (associated(ptr))
				if (path == ptr%id) then
					ok=.true.
					exit
				else
					if (present(ptr2)) ptr2=>ptr
					ptr=>next_section(ptr)
				end if
			end do
		end if
	end function

	recursive subroutine read_section(sect)
		type(section_t) :: sect
		logical :: sarg

        character(6) :: str, str2
        integer(SP) :: nkw, nsec, i, iarg
		type(keyword_t), pointer :: kw
		type(section_t), pointer :: next

		read(lun,*) str, sect%id, sect%nsect, sect%set
		if (sane('SECT', str)) stop

		read(lun,*) str, sarg, str2, sect%nkw
		if (sane('ARG', str)) stop
		if (sane('KW', str2)) stop
		
		if (sarg) then
			call new_keyword(sect%arg)
			call read_keyword(sect%arg)
		end if

		if (sect%nkw > 0) then
			call new_keyword(sect%kw)
			call read_keyword(sect%kw)
			kw=>sect%kw
			do i=2, sect%nkw
				call new_keyword(kw%next)
				call read_keyword(kw%next)
				kw=>kw%next
			end do
		end if

		if (sect%nsect > 0) then
			call new_section(sect%sects)
			call read_section(sect%sects)
			next=>sect%sects
			do i=2, sect%nsect
				call new_section(next%next)
				call read_section(next%next)
				next=>next%next
			end do
		end if
	end subroutine

	subroutine read_keyword(kw) 
		type(keyword_t) :: kw

        character(MAXID) :: typ
        integer(SP) :: n

		read(lun,*) typ, kw%id, n, kw%set
		call new_kword(kw%key)
		select case(typ)
			case('INT')
				kw%typ=KW_INT
				if (n > 0) then
!                    call read_int_kw(kw%key,n)
					read(lun,*) kw%key%ival
					kw%key%set=.true.
				end if
			case('INT_ARRAY')
				kw%typ=KW_IVEC
				if (n > 0) then
					call read_int_kw(kw%key,n)
					kw%key%set=.true.
				end if
			case('DBL')
				kw%typ=KW_DBL
				if (n > 0) then
!                    call read_dbl_kw(kw%key,n)
					read(lun,*) kw%key%dval
					kw%key%set=.true.
				end if
			case('DBL_ARRAY')
				kw%typ=KW_DVEC
				if (n > 0) then
					call read_dbl_kw(kw%key,n)
					kw%key%set=.true.
				end if
			case('BOOL')
				kw%typ=KW_BOOL
				if (n > 0) then 
!                    call read_bool_kw(kw%key,n)
					read(lun,*) kw%key%bool
					kw%key%set=.true.
				end if
			case('BOOL_ARRAY')
				kw%typ=KW_LVEC
				if (n > 0) then 
					call read_bool_kw(kw%key,n)
					kw%key%set=.true.
				end if
			case('STR')
				kw%typ=KW_STR
				if (n == -1) then ! empty string
					allocate(kw%key%str(1))
					memtrack=memtrack+1
					kw%key%str(1)=''
					kw%key%set=.true.
				else if (n > 0) then
					call read_str_kw(kw%key,n)
					kw%key%set=.true.
				else
					nullify(kw%key%str)
				end if
			case default
				call msg_error('invalid type: ' // trim(kw%id) // ' -> ' // typ)
				stop
		end select
!        call print_kw(kw)
	end subroutine

	subroutine read_int_kw(kw,n)
		type(kword_t) :: kw
		integer(SP) :: n

!        if (n > 1) then
			allocate(kw%ivec(n))
			memtrack=memtrack+1
			read(lun,*) kw%ivec
!        else
!            read(lun,*) kw%ival
!        end if
		
	end subroutine

	subroutine read_dbl_kw(kw,n)
		type(kword_t) :: kw
		integer(SP) :: n

!        if (n > 1) then
			allocate(kw%dvec(n))
			memtrack=memtrack+1
			read(lun,*) kw%dvec
!        else
!            read(lun,*) kw%dval
!        end if
	end subroutine

	subroutine read_str_kw(kw,n)
		type(kword_t) :: kw
		integer(SP) :: n

		integer(SP) :: strmx

!        read(lun,'(a)') strmx
		allocate(kw%str(n))
		memtrack=memtrack+1
		read(lun,'(a)') kw%str
	end subroutine

	subroutine read_bool_kw(kw,n)
		type(kword_t) :: kw
		integer(SP) :: n

!        if (n > 1) then
			allocate(kw%lvec(n))
			memtrack=memtrack+1
			read(lun,*) kw%lvec
!        else
			read(lun,*) kw%bool
!        end if
	end subroutine

	function sane(str1, str2) result(err)
		character(*), intent(in) :: str1, str2
		logical :: err

		err=.false.
		if (str1 /= str2) then
			call msg_error('parse error: ' //trim(str2)// '!=' //trim(str1))
			err=.true.
		end if
	end function

	recursive subroutine print_section(sect)
		type(section_t) :: sect

		type(keyword_t), pointer :: kw
		type(section_t), pointer :: ss

		print *, 'SECTION: ', sect%id, sect%nsect, sect%nkw
		if (associated(sect%arg)) then
			print *, 'ARG:'
			call print_kw(sect%arg)
		end if

		if (associated(sect%kw)) then
			kw=>sect%kw
			do while (associated(kw))
				call print_kw(kw)
				kw=>kw%next
			end do
		end if

		if (associated(sect%sects)) then
			call print_section(sect%sects)
		end if

		if (associated(sect%next)) then
			call print_section(sect%next)
		end if
	end subroutine

    subroutine print_kw(kw)
        type(keyword_t) :: kw

		integer :: i

		if (.not.kw%set) then
			print *, trim(kw%id), 'not set.'
			return
		end if
		select case (kw%typ)
			case(KW_INT)
				print *, trim(kw%id), kw%key%ival
			case(KW_IVEC)
				print *, trim(kw%id), kw%key%ivec
			case(KW_DBL)
				print *, trim(kw%id), kw%key%dval
			case(KW_DVEC)
				print *, trim(kw%id), kw%key%dvec
			case(KW_BOOL)
				print *, trim(kw%id), kw%key%bool
			case(KW_LVEC)
				print *, trim(kw%id), kw%key%lvec
			case(KW_STR)
				print *, trim(kw%id), ':'
				do i=1,size(kw%key%str)
					print *, trim(kw%key%str(i))
				end do
			case(KW_NONE)
				call msg_warn('key undefined')
			case default
				call msg_warn('dang! no such key type' // xstr(kw%typ))
		end select
	end subroutine

	subroutine push_section(self, sect)
		type(getkw_t), target :: self
		character(*) :: sect

		type(section_t), pointer :: ptr
		logical :: ok

		if (self%stidx < MAX_STACK) then
			self%stidx=self%stidx+1
			self%stack(self%stidx)%sect=>self%active
			if (sect == '.') then
				self%active=>self%stack(1)%sect
				return
			end if
			ok=find_sect(self%active, sect, ptr)
			if ( .not.ok ) then
				call msg_error('Invalid section: ' // sect)
				stop
			end if
			self%active=>ptr
		else if (strict) then 
			call msg_error('push_section: stack overflow!')
			stop
		else
			call msg_warn('push_section: stack overflow!')
		end if
			
	end subroutine 

	subroutine pop_section(self)
		type(getkw_t), target :: self
		integer(4) :: error

		if (self%stidx > 0) then
			self%active=>self%stack(self%stidx)%sect
			self%stidx=self%stidx-1
		else if (strict) then
			call msg_error('pop_section: stack underflow!')
			stop
		else
			call msg_warn('pop_section: stack underflow!')
		end if
	end subroutine 

	subroutine reset_active_section(self)
		type(getkw_t), target :: self
		
		self%active=>self%main
	end subroutine 

	subroutine validate(error, str)
		integer(4), intent(in) :: error
		character(*), optional :: str

		if (error == 0) then
			if (verbose) then
				if (present(str)) then
					print *, '<<< Keyword lookup failed: ', str, ' >>>'
				else
					print *, '<<< Keyword lookup failed! >>>'
				end if
			end if
			if ( strict ) stop
		end if
	end subroutine 

	function keyword_is_set(self, key) result(ok)
		type(getkw_t) :: self
		character(*), intent(in) :: key
		logical :: ok

		type(keyword_t), pointer :: ptr

		ok=findkw(self%active, key, ptr)
		if (.not.ok) then
			call msg_warn('getkw: no such key: ' // trim(key))
		else
			ok=ptr%set
		end if
	end function 

	function section_is_set(self, key) result(ok)
		type(getkw_t) :: self
		character(*), intent(in) :: key
		logical :: ok

		type(section_t), pointer :: ptr

		ok=find_sect(self%active, key, ptr)
		if (.not.ok) then
			call msg_warn('no such section: ' // trim(key))
		else
			ok=ptr%set
		end if
	end function 

	function has_keyword(self, key) result(ok)
		type(getkw_t) :: self
		character(*), intent(in) :: key
		logical :: ok

		type(keyword_t), pointer :: ptr

		ok=findkw(self%active, key, ptr)
	end function 

	function next_keyword(key) result( new)
		type(keyword_t) :: key
		type(keyword_t), pointer :: new

		new=>key%next
	end function

	function next_kword(key) result( new)
		type(kword_t) :: key
		type(kword_t), pointer :: new

		new=>key%next
	end function

	function next_section(sect) result( new)
		type(section_t) :: sect
		type(section_t), pointer :: new

		new=>sect%next
	end function

	function last_keyword(key) result( new)
		type(keyword_t), pointer :: key
		type(keyword_t), pointer :: new

		new=>key
		do while (associated(new%next))
			new=>new%next
		end do
	end function

	function last_kword(key) result( new)
		type(kword_t), pointer :: key
		type(kword_t), pointer :: new

!        new=>key%next
		new=>key
!        if (.not.associated(key)) return
		do while (associated(new%next))
			new=>new%next
		end do
	end function

	function last_section(sect) result( new)
		type(section_t), pointer :: sect
		type(section_t), pointer :: new

		new=>sect
!        if (.not.associated(sect)) return
		do while (associated(new%next))
			new=>new%next
		end do
	end function

	subroutine set_verbose(v)
		logical, intent(in) :: v
		verbose=v
	end subroutine

	subroutine set_strict(s)
		logical, intent(in) :: s
		strict=s
	end subroutine

	subroutine print_tree(t)
		type(getkw_t) :: t

		type(section_t), pointer :: s

		s=>t%main
		call print_section(s)

	end subroutine

end module
	
