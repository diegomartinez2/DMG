!
! Simple readline functionality
!
subroutine showline(histline, hpos, history)
  use curses
  integer, parameter :: MAXLINES = 20
  integer, intent(in) :: histline
  integer, intent(in) :: hpos
  character (len=4) :: curline
  character (len=20) :: xpos, ypos
  character (len=80), dimension(MAXLINES), intent(in) :: history
  integer :: ierr
!
  write(curline,'(i4)') histline
! ierr=erase()
  ierr=move(1,1)
  ierr=clrtoeol()
  ierr=mvaddstr(1, 1, trim(curline) // '> ' // C_NULL_CHAR)
  ierr=addstr(history(histline)(1:hpos-1) // C_NULL_CHAR)
  ierr=attrset(curses_a_reverse)
! ierr=attrset(262144)
  ierr=addch(ichar(history(histline)(hpos:hpos)))
  ierr=attrset(curses_a_normal)
  ierr=addstr(history(histline)(hpos+1:len_trim(history(histline))) // C_NULL_CHAR)
  ierr=refresh()
end subroutine showline

program readline
  use curses
  integer, parameter :: ascii_bs = 8, ascii_lf = 10, ascii_cr = 13
  integer, parameter :: MAXLINES = 20
  character (len=80), dimension(MAXLINES) :: history
  integer :: histline, hpos, linlen, mode, nhist
  integer :: ierr
  type (C_PTR) :: iwin

  istat=0
  mode=1
  iwin = C_NULL_PTR
  iwin = initscr()
  histline=1
  history(histline)=' '
  nhist=1
  hpos=1
  linlen=0
  call showline(histline, hpos, history)
  do
    ierr=noecho()
    ierr=raw()
    ierr=curs_set(0)
    ierr=keypad(stdscr, curses_true)
    key=getch()
!   if (key > 0) write(*,*) 'key= ', key
    if (key == ascii_lf .or. key == ascii_cr .or. key==curses_key_enter) then
      if (history(histline) == 'quit') exit
      histline=histline+1
      if (histline > MAXLINES) histline=MAXLINES
      history(histline)=' '
      hpos=1
      nhist=nhist+1
      call showline(histline, hpos, history)
    else if (key == curses_key_up) then
      histline=histline-1
      if (histline == 0) histline=1
      hpos=1
      call showline(histline, hpos, history)
    else if (key == curses_key_left) then
      hpos=hpos-1
      if (hpos == 0) hpos=1
      call showline(histline, hpos, history)
    else if (key == curses_key_right) then
      hpos=hpos+1
      if (hpos > 80) hpos=80
      call showline(histline, hpos, history)
    else if (key == curses_key_down) then
      histline=histline+1
      if (histline > nhist) histline=nhist
      hpos=1
      call showline(histline, hpos, history)
    else if (key == 14) then
      linlen=len_trim(history(histline))
      do
        if (hpos > linlen) then
          hpos=linlen
          exit
        end if
        if (history(histline)(hpos:hpos) == ' ') exit
        hpos=hpos+1
      end do
      do
        if (hpos > linlen) then
          hpos=linlen
          exit
        end if
        if (history(histline)(hpos:hpos) /= ' ') then
          exit
        end if
        hpos=hpos+1
      end do
      call showline(histline, hpos, history)
    else if (key == 16) then
      hpos=hpos-1
      do
        if (hpos < 1) then
          hpos=1
          exit
        end if
        if (history(histline)(hpos:hpos) /= ' ') exit
        hpos=hpos-1
      end do
      do
        if (hpos < 1) then
          hpos=1
          exit
        end if
        if (history(histline)(hpos:hpos) == ' ') then
          hpos=min(hpos+1, len_trim(history(histline)))
          exit
        end if
        hpos=hpos-1
      end do
      call showline(histline, hpos, history)
    else if (key == curses_key_ic) then
      mode=3-mode
    else if (key == curses_key_dc .or. key == ascii_bs) then
      if (key == ascii_bs) then
        hpos=hpos-1
        if (hpos == 0) hpos=1
      end if
      linlen=len_trim(history(histline))
      history(histline)(hpos:linlen)=history(histline)(hpos+1:linlen) // ' '
      call showline(histline, hpos, history)
    else if (key > 31 .and. key < 256) then
      if (mode == 1) then
        linlen=len_trim(history(histline))
        history(histline)(hpos+1:linlen+1)=history(histline)(hpos:linlen)
      end if
      history(histline)(hpos:hpos)=achar(key)
      hpos=hpos+1
      if (hpos > 80) hpos=80
      call showline(histline, hpos, history)
    else if (key == 3) then
      exit
    end if
  end do
  ierr = endwin()
end program readline
