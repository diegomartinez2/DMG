module ncurses
  use iso_c_binding
  implicit none

  integer (C_SIGNED_CHAR), parameter :: curses_true = 1
  integer (C_SIGNED_CHAR), parameter :: curses_false = 0
  integer, parameter :: curses_err = -1
  integer, parameter :: curses_ok  =  0
! ---------------------------------------------------------------------
!
!  PDCurses Mouse Interface -- SYSVR4, with extensions
!
  type mouse_data
    integer(C_INT) :: x        ! absolute column, 0 based, measured in characters
    integer(C_INT) :: y        ! absolute row, 0 based, measured in characters
    integer(C_SHORT), dimension(3) :: button ! state of each button
    integer(C_INT) :: changes  ! flags indicating what has changed with the mouse
  end type mouse_data

  type (mouse_data) :: mouse_status

!------------------------------------------------------------------------
!
!   PDCurses Structure Definitions
!
! definition of a window
! type window_data
!   integer (C_INT) :: cury, curx    ! current pseudo-cursor
!   integer (C_INT) :: maxy, maxx    ! max window coordinates
!   integer (C_INT) :: begy, begx    ! origin on screen
!   integer (C_INT) :: flags         ! window properties
!   integer (C_LONG) ::  attrs       ! standard attributes and colors
!   integer (C_LONG) ::  bkgd        ! background, normally blank
!   integer (C_SIGNED_CHAR) :: clear      ! causes clear at next refresh
!   integer (C_SIGNED_CHAR) :: leaveit    ! leaves cursor where it is
!   integer (C_SIGNED_CHAR) :: scroll     ! allows window scrolling
!   integer (C_SIGNED_CHAR) :: nodelay    ! input character wait flag
!   integer (C_SIGNED_CHAR) :: immed      ! immediate update flag
!   integer (C_SIGNED_CHAR) :: sync       ! synchronise window ancestors
!   integer (C_SIGNED_CHAR) :: use_keypad ! flags keypad key mode active
!   type (C_PTR) :: y                ! pointer to line pointer array
!   type (C_PTR) :: firstch          ! first changed character in line
!   type (C_PTR) :: lastch           ! last changed character in line
!   integer (C_INT) :: tmarg         ! top of scrolling region
!   integer (C_INT) :: bmarg         ! bottom of scrolling region
!   integer (C_INT) :: delayms       ! milliseconds of delay for getch()
!   integer (C_INT) :: parx, pary    ! coords relative to parent (0,0)
!   type (C_PTR) :: parent           ! subwin pointer to parent win
! end type window_data
!
  integer, parameter :: curses_button_released         =  0
  integer, parameter :: curses_button_pressed          =  1
  integer, parameter :: curses_button_clicked          =  2
  integer, parameter :: curses_button_double_clicked   =  3
  integer, parameter :: curses_button_triple_clicked   =  4
  integer, parameter :: curses_button_moved            =  5
  integer, parameter :: curses_wheel_scrolled          =  6
  integer, parameter :: curses_button_action_mask      =  7
  integer, parameter :: curses_pdc_button_shift        =  8
  integer, parameter :: curses_pdc_button_control      = 16
  integer, parameter :: curses_pdc_button_alt          = 32
  integer, parameter :: curses_button_modifier_mask    = 56

  integer, parameter :: curses_color_black     = 0
  integer, parameter :: curses_color_red       = 1
  integer, parameter :: curses_color_green     = 2
  integer, parameter :: curses_color_yellow    = 3
  integer, parameter :: curses_color_blue      = 4
  integer, parameter :: curses_color_magenta   = 5
  integer, parameter :: curses_color_cyan      = 6
  integer, parameter :: curses_color_white     = 7

  integer, parameter :: curses_key_code_yes  = 256 !  If get_wch() gives a key code

  integer, parameter :: curses_key_break     = 257 ! Not on PC KBD
  integer, parameter :: curses_key_down      = 258 ! Down arrow key
  integer, parameter :: curses_key_up        = 259 ! Up arrow key
  integer, parameter :: curses_key_left      = 260 ! Left arrow key
  integer, parameter :: curses_key_right     = 261 ! Right arrow key
  integer, parameter :: curses_key_home      = 262 ! home key
  integer, parameter :: curses_key_backspace = 263 ! not on pc
  integer, parameter :: curses_key_f0        = 264 ! function keys; 64 reserved

  integer, parameter :: curses_key_dl        = 328 ! delete line
  integer, parameter :: curses_key_il        = 329 ! insert line
  integer, parameter :: curses_key_dc        = 330 ! delete character
  integer, parameter :: curses_key_ic        = 331 ! insert char or enter ins mode
  integer, parameter :: curses_key_eic       = 332 ! exit insert char mode
  integer, parameter :: curses_key_clear     = 333 ! clear screen
  integer, parameter :: curses_key_eos       = 334 ! clear to end of screen
  integer, parameter :: curses_key_eol       = 335 ! clear to end of line
  integer, parameter :: curses_key_sf        = 336 ! scroll 1 line forward
  integer, parameter :: curses_key_sr        = 337 ! scroll 1 line back (reverse)
  integer, parameter :: curses_key_npage     = 338 ! next page
  integer, parameter :: curses_key_ppage     = 339 ! previous page
  integer, parameter :: curses_key_stab      = 340 ! set tab
  integer, parameter :: curses_key_ctab      = 341 ! clear tab
  integer, parameter :: curses_key_catab     = 342 ! clear all tabs
  integer, parameter :: curses_key_enter     = 343 ! enter or send (unreliable)
  integer, parameter :: curses_key_sreset    = 344 ! soft/reset (partial/unreliable)
  integer, parameter :: curses_key_reset     = 345 ! reset/hard reset (unreliable)
  integer, parameter :: curses_key_print     = 346 ! print/copy
  integer, parameter :: curses_key_ll        = 347 ! home down/bottom (lower left)
  integer, parameter :: curses_key_a1        = 348 ! a1
  integer, parameter :: curses_key_a3        = 349 ! a3
  integer, parameter :: curses_key_b2        = 350 ! b2
  integer, parameter :: curses_key_c1        = 351 ! c1
  integer, parameter :: curses_key_c3        = 352 ! c3
  integer, parameter :: curses_key_btab      = 353 ! Back tab key
  integer, parameter :: curses_key_beg       = 354 ! beg(inning) key
  integer, parameter :: curses_key_cancel    = 355 ! cancel key
  integer, parameter :: curses_key_close     = 356 ! close key
  integer, parameter :: curses_key_command   = 357 ! cmd (command) key
  integer, parameter :: curses_key_copy      = 358 ! copy key
  integer, parameter :: curses_key_create    = 359 ! create key
  integer, parameter :: curses_key_end       = 360 ! end key
  integer, parameter :: curses_key_exit      = 361 ! exit key
  integer, parameter :: curses_key_find      = 362 ! find key
  integer, parameter :: curses_key_help      = 363 ! help key
  integer, parameter :: curses_key_mark      = 364 ! mark key
  integer, parameter :: curses_key_message   = 365 ! message key
  integer, parameter :: curses_key_move      = 366 ! move key
  integer, parameter :: curses_key_next      = 367 ! next object key
  integer, parameter :: curses_key_open      = 368 ! open key
  integer, parameter :: curses_key_options   = 369 ! options key
  integer, parameter :: curses_key_previous  = 370 ! previous object key
  integer, parameter :: curses_key_redo      = 371 ! redo key
  integer, parameter :: curses_key_reference = 372 ! ref(erence) key
  integer, parameter :: curses_key_refresh   = 373 ! refresh key
  integer, parameter :: curses_key_replace   = 374 ! replace key
  integer, parameter :: curses_key_restart   = 375 ! restart key
  integer, parameter :: curses_key_resume    = 376 ! resume key
  integer, parameter :: curses_key_save      = 377 ! save key
  integer, parameter :: curses_key_sbeg      = 378 ! shifted beginning key
  integer, parameter :: curses_key_scancel   = 379 ! shifted cancel key
  integer, parameter :: curses_key_scommand  = 380 ! shifted command key
  integer, parameter :: curses_key_scopy     = 381 ! shifted copy key
  integer, parameter :: curses_key_screate   = 382 ! shifted create key
  integer, parameter :: curses_key_sdc       = 383 ! shifted delete char key
  integer, parameter :: curses_key_sdl       = 384 ! shifted delete line key
  integer, parameter :: curses_key_select    = 385 ! select key
  integer, parameter :: curses_key_send      = 386 ! shifted end key
  integer, parameter :: curses_key_seol      = 387 ! shifted clear line key
  integer, parameter :: curses_key_sexit     = 388 ! shifted exit key
  integer, parameter :: curses_key_sfind     = 389 ! shifted find key
  integer, parameter :: curses_key_shelp     = 390 ! shifted help key
  integer, parameter :: curses_key_shome     = 391 ! shifted home key
  integer, parameter :: curses_key_sic       = 392 ! shifted input key
  integer, parameter :: curses_key_sleft     = 393 ! shifted left arrow key
  integer, parameter :: curses_key_smessage  = 394 ! shifted message key
  integer, parameter :: curses_key_smove     = 395 ! shifted move key
  integer, parameter :: curses_key_snext     = 396 ! shifted next key
  integer, parameter :: curses_key_soptions  = 397 ! shifted options key
  integer, parameter :: curses_key_sprevious = 398 ! shifted prev key
  integer, parameter :: curses_key_sprint    = 399 ! shifted print key
  integer, parameter :: curses_key_sredo     = 400 ! shifted redo key
  integer, parameter :: curses_key_sreplace  = 401 ! shifted replace key
  integer, parameter :: curses_key_sright    = 402 ! shifted right arrow
  integer, parameter :: curses_key_srsume    = 403 ! shifted resume key
  integer, parameter :: curses_key_ssave     = 404 ! shifted save key
  integer, parameter :: curses_key_ssuspend  = 405 ! shifted suspend key
  integer, parameter :: curses_key_sundo     = 406 ! shifted undo key
  integer, parameter :: curses_key_suspend   = 407 ! suspend key
  integer, parameter :: curses_key_undo      = 408 ! undo key
  integer, parameter :: curses_key_mouse     = 409 ! "mouse" key
  integer, parameter :: curses_key_resize    = 410 ! window resize
  integer, parameter :: curses_key_event     = 411 ! event key
  integer, parameter :: curses_key_max       = 511 ! undo key
!
! The available non-color attributes are bold, underline, invisible,
! right-line, left-line, protect, reverse and blink. 256 color pairs (8
! bits), 8 bits for other attributes, and 16 bits for character data.
!
  integer(C_LONG), parameter :: curses_a_normal     = 0
  integer(C_LONG), parameter :: curses_a_altcharset = 65536
  integer(C_LONG), parameter :: curses_a_rightline  = 268435456
  integer(C_LONG), parameter :: curses_a_leftline   = 67108864
  integer(C_LONG), parameter :: curses_a_invis      = 8388608
  integer(C_LONG), parameter :: curses_a_underline  = 131072
  integer(C_LONG), parameter :: curses_a_reverse    = 262144
  integer(C_LONG), parameter :: curses_a_blink      = 524288
  integer(C_LONG), parameter :: curses_a_bold       = 2097152

  integer(C_LONG), parameter :: curses_a_attributes = -256
  integer(C_LONG), parameter :: curses_a_chartext   = 255
  integer(C_LONG), parameter :: curses_a_color      = 65280

  integer(C_LONG), parameter :: curses_a_italic     = 65536
  integer(C_LONG), parameter :: curses_a_protect    = 16777216

!
!  VT100-compatible symbols -- box chars
!
  integer(C_LONG), parameter :: acs_ulcorner = 65644
  integer(C_LONG), parameter :: acs_llcorner = 65645
  integer(C_LONG), parameter :: acs_urcorner = 65643
  integer(C_LONG), parameter :: acs_lrcorner = 65642
  integer(C_LONG), parameter :: acs_rtee     = 65653
  integer(C_LONG), parameter :: acs_ltee     = 65652
  integer(C_LONG), parameter :: acs_btee     = 65654
  integer(C_LONG), parameter :: acs_ttee     = 65655
  integer(C_LONG), parameter :: acs_hline    = 65649
  integer(C_LONG), parameter :: acs_vline    = 65656
  integer(C_LONG), parameter :: acs_plus     = 65646
!
! VT100-compatible symbols -- other
!
  integer(C_LONG), parameter :: acs_s1       = 65647
  integer(C_LONG), parameter :: acs_s9       = 65651
  integer(C_LONG), parameter :: acs_diamond  = 65632
  integer(C_LONG), parameter :: acs_ckboard  = 65633
  integer(C_LONG), parameter :: acs_degree   = 65638
  integer(C_LONG), parameter :: acs_plminus  = 65639
  integer(C_LONG), parameter :: acs_bullet   = 65662
!
! Teletype 5410v1 symbols -- these are defined in SysV curses, but
! are not well-supported by most terminals. Stick to VT100 characters
! for optimum portability.
!
  integer(C_LONG), parameter :: acs_larrow   = 65580
  integer(C_LONG), parameter :: acs_rarrow   = 65579
  integer(C_LONG), parameter :: acs_darrow   = 65582
  integer(C_LONG), parameter :: acs_uarrow   = 65581
  integer(C_LONG), parameter :: acs_board    = 65640
  integer(C_LONG), parameter :: acs_lantern  = 65641
  integer(C_LONG), parameter :: acs_block    = 65584
!----------------------------------------------------------
!
! PDCurses External Variables
!
  integer (C_INT), bind(C, name='LINES') :: LINES
  integer (C_INT), bind(C, name='COLS') :: COLS
! windows
  type    (C_PTR), bind(C, name='stdscr') :: stdscr
  type    (C_PTR), bind(C, name='curscr') :: curscr
! screen
  type    (C_PTR), bind(C, name='sp') :: sp
  integer (C_INT), bind(C, name='COLORS') :: COLORS
  integer (C_INT), bind(C, name='COLOR_PAIRS') :: COLOR_PAIRS
  integer (C_INT), bind(C, name='TABSIZE') :: TABSIZE
! integer(C_LONG) :: acs_map
!
!----------------------------------------------------------------------
!
!   PDCurses Function Declarations
!
! Standard
!
! addch, addchnstr, addchstr, addnstr, addstr,
! attroff, attron, attrset, attr_get, attr_off,
! attr_on, attr_set,
! baudrate, beep, bkgd, bkgdset, border, box,
! can_change_color, cbreak, chgat, clearok, clear, clrtobot,
! clrtoeol, color_content, color_set, copywin, curs_set,
! def_prog_mode, def_shell_mode, delay_output,
! delch, deleteln, delscreen, delwin,
! derwin, doupdate, dupwin,
! echochar, echo, endwin, erasechar, erase,
! filter, flash, flushinp,
! getbkgd, getnstr, getstr, getwin,
! halfdelay, has_colors, has_ic, has_il, hline,
! idcok, idlok, immedok, inchnstr, inchstr, inch,
! init_color, init_pair, initscr, innstr,
! insch, insdelln, insertln, insnstr, insstr, instr,
! intrflush, isendwin, is_linetouched, is_wintouched,
! keyname, keypad, killchar, leaveok, longname,
! meta, move, mvaddch, mvaddchnstr, mvaddchstr,
! mvaddnstr, mvaddstr, mvchgat, mvcur, mvdelch,
! mvderwin, mvgetch, mvgetnstr, mvgetstr, mvhline, mvinch,
! mvinchnstr, mvinchstr, mvinnstr, mvinsch, mvinsnstr, mvinsstr,
! mvinstr, mvprintw, mvscanw, mvvline, mvwaddchnstr,
! mvwaddchstr, mvwaddch, mvwaddnstr, mvwaddstr,
! mvwchgat, mvwdelch, mvwgetch, mvwgetnstr, mvwgetstr,
! mvwhline, mvwinchnstr, mvwinchstr, mvwinch, mvwinnstr,
! mvwinsch, mvwinsnstr, mvwinsstr, mvwinstr, mvwin,
! mvwprintw, mvwscanw, mvwvline,
! napms, newpad, newterm, newwin, nl,
! nocbreak, nodelay, noecho, nonl, noqiflush, noraw, notimeout,
! overlay, overwrite,
! pair_content, pechochar, pnoutrefresh, prefresh,
! printw, putwin,
! qiflush, raw, redrawwin, refresh,
! reset_prog_mode, reset_shell_mode, resetty, ripoffline,
! savetty, scanw, scr_dump, scr_init, scr_restore, scr_set,
! scrl, scroll, scrollok, set_term, setscrreg,
! slk_attroff, slk_attr_off, slk_attron, slk_attr_on, slk_attrset,
! slk_attr_set, slk_clear, slk_color, slk_init, slk_label,
! slk_noutrefresh, slk_refresh, slk_restore, slk_set, slk_touch,
! standend, standout, start_color,
! subpad, subwin, syncok,
! termattrs, term_attrs, termname, timeout,
! touchline, touchwin, typeahead,
! untouchwin, use_env,
! vidattr, vid_attr, vidputs, vid_puts,
! vline, vw_printw, vwprintw, vw_scanw, vwscanw,
! waddchnstr, waddchstr, waddch, waddnstr, waddstr,
! wattroff, wattron, wattrset, wattr_get, wattr_off, wattr_on,
! wattr_set, wbkgdset, wbkgd, wborder, wchgat, wclear,
! wclrtobot, wclrtoeol, wcolor_set, wcursyncup, wdelch,
! wdeleteln, wechochar, werase, wgetch, wgetnstr,
! wgetstr, whline, winchnstr, winchstr, winch,
! winnstr, winsch, winsdelln, winsertln, winsnstr, winsstr,
! winstr, wmove, wnoutrefresh, wprintw, wredrawln, wrefresh,
! wscanw, wscrl, wsetscrreg, wstandend, wstandout, wsyncdown,
! wsyncup, wtimeout, wtouchln, wvline
!
  interface
!         int addch(const chtype ch);
    function addch(ch) result (ires) bind(C, name='addch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: ch
    end function addch
!         int waddch(WINDOW *win, const chtype ch);
    function waddch(win,ch) result (ires) bind(C, name='waddch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value, intent(in) :: ch
    end function waddch
!         int mvaddch(int y, int x, const chtype ch);
    function mvaddch(y,x,ch) result (ires) bind(C, name='mvaddch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value, intent(in) :: ch
    end function mvaddch
!         int mvwaddch(WINDOW *win, int y, int x, const chtype ch);
    function mvwaddch(win,y,x,ch) result (ires) bind(C, name='mvwaddch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value, intent(in) :: ch
    end function mvwaddch
!         int echochar(const chtype ch);
    function echochar(ch) result (ires) bind(C, name='echochar')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: ch
    end function echochar
!         int wechochar(WINDOW *win, const chtype ch);
    function wechochar(win,ch) result (ires) bind(C, name='wechochar')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value, intent(in) :: ch
    end function wechochar

!         int addrawch(chtype ch);
    function addrawch(ch) result (ires) bind(C, name='addrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function addrawch
!         int waddrawch(WINDOW *win, chtype ch);
    function waddrawch(win,ch) result (ires) bind(C, name='waddrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function waddrawch
!         int mvaddrawch(int y, int x, chtype ch);
    function mvaddrawch(y,x,ch) result (ires) bind(C, name='mvaddrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvaddrawch
!         int mvwaddrawch(WINDOW *win, int y, int x, chtype ch);
    function mvwaddrawch(win,y,x,ch) result (ires) bind(C, name='mvwaddrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvwaddrawch

!         int addchstr(const chtype *ch);
    function addchstr(ch) result (ires) bind(C, name='addchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), intent(in), value :: ch
    end function addchstr
!         int addchnstr(const chtype *ch, int n);
    function addchnstr(ch,n) result (ires) bind(C, name='addchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function addchnstr
!         int waddchstr(WINDOW *win, const chtype *ch);
    function waddchstr(win,ch) result (ires) bind(C, name='waddchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), intent(in), value :: ch
    end function waddchstr
!         int waddchnstr(WINDOW *win, const chtype *ch, int n);
    function waddchnstr(win,ch,n) result (ires) bind(C, name='waddchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function waddchnstr
!         int mvaddchstr(int y, int x, const chtype *ch);
    function mvaddchstr(y,x,ch) result (ires) bind(C, name='mvaddchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
    end function mvaddchstr
!         int mvaddchnstr(int y, int x, const chtype *ch, int n);
    function mvaddchnstr(y,x,ch,n) result (ires) bind(C, name='mvaddchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function mvaddchnstr
!         int mvwaddchstr(WINDOW *, int y, int x, const chtype *ch);
    function mvwaddchstr(win,y,x,ch) result (ires) bind(C, name='mvwaddchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
    end function mvwaddchstr
!         int mvwaddchnstr(WINDOW *, int y, int x, const chtype *ch, int n);
    function mvwaddchnstr(win,y,x,ch,n) result (ires) bind(C, name='mvwaddchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function mvwaddchnstr

!         int addstr(const char *str);
    function addstr(str) result (ires) bind(C, name='addstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
    end function addstr
!         int addnstr(const char *str, int n);
    function addnstr(str,n) result (ires) bind(C, name='addnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function addnstr
!         int waddstr(WINDOW *win, const char *str);
    function waddstr(win,str) result (ires) bind(C, name='waddstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
    end function waddstr
!         int waddnstr(WINDOW *win, const char *str, int n);
    function waddnstr(win,str,n) result (ires) bind(C, name='waddnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function waddnstr
!         int mvaddstr(int y, int x, const char *str);
    function mvaddstr(y,x,str) result (ires) bind(C, name='mvaddstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvaddstr
!         int mvaddnstr(int y, int x, const char *str, int n);
    function mvaddnstr(y,x,str,n) result (ires) bind(C, name='mvaddnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvaddnstr
!         int mvwaddstr(WINDOW *win, int y, int x, const char *str);
    function mvwaddstr(win,y,x,str) result (ires) bind(C, name='mvwaddstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvwaddstr
!         int mvwaddnstr(WINDOW *win, int y, int x, const char *str, int n);
    function mvwaddnstr(win,y,x,str,n) result (ires) bind(C, name='mvwaddnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvwaddnstr

!         int attroff(chtype attrs);
    function attroff(attrs) result (ires) bind(C, name='attroff')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
    end function attroff
!         int wattroff(WINDOW *win, chtype attrs);
    function wattroff(win,attrs) result (ires) bind(C, name='wattroff')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
    end function wattroff
!         int attron(chtype attrs);
    function attron(attrs) result (ires) bind(C, name='attron')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
    end function attron
!         int wattron(WINDOW *win, chtype attrs);
    function wattron(win,attrs) result (ires) bind(C, name='wattron')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
    end function wattron
!         int attrset(chtype attrs);
    function attrset(attrs) result (ires) bind(C, name='attrset')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
    end function attrset
!         int wattrset(WINDOW *win, chtype attrs);
    function wattrset(win,attrs) result (ires) bind(C, name='wattrset')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
    end function wattrset
!         int standend(void);
    function standend() result (ires) bind(C, name='standend')
      use iso_c_binding
      integer (C_INT) :: ires
    end function standend
!         int wstandend(WINDOW *win);
    function wstandend(win) result (ires) bind(C, name='wstandend')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wstandend
!         int standout(void);
    function standout() result (ires) bind(C, name='standout')
      use iso_c_binding
      integer (C_INT) :: ires
    end function standout
!         int wstandout(WINDOW *win);
    function wstandout(win) result (ires) bind(C, name='wstandout')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wstandout

!         int color_set(short color_pair, void *opts);
    function color_set(color_pair,opts) result (ires) bind(C, name='color_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function color_set
!         int wcolor_set(WINDOW *win, short color_pair, void *opts);
    function wcolor_set(win,color_pair,opts) result (ires) bind(C, name='wcolor_set')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function wcolor_set

!         int attr_get(attr_t *attrs, short *color_pair, void *opts);
    function attr_get(attrs,color_pair,opts) result (ires) bind(C, name='attr_get')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: attrs
      type (C_PTR), value :: color_pair
      type (C_PTR), value :: opts
    end function attr_get
!         int attr_off(attr_t attrs, void *opts);
    function attr_off(attrs,opts) result (ires) bind(C, name='attr_off')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function attr_off
!         int attr_on(attr_t attrs, void *opts);
    function attr_on(attrs,opts) result (ires) bind(C, name='attr_on')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function attr_on
!         int attr_set(attr_t attrs, short color_pair, void *opts);
    function attr_set(attrs,color_pair,opts) result (ires) bind(C, name='attr_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function attr_set
!         int wattr_get(WINDOW *win, attr_t *attrs, short *color_pair, void *opts);
    function wattr_get(win,attrs,color_pair,opts) result (ires) bind(C, name='wattr_get')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), value :: attrs
      type (C_PTR), value :: color_pair
      type (C_PTR), value :: opts
    end function wattr_get
!         int wattr_off(WINDOW *win, attr_t attrs, void *opts);
    function wattr_off(win,attrs,opts) result (ires) bind(C, name='wattr_off')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function wattr_off
!         int wattr_on(WINDOW *win, attr_t attrs, void *opts);
    function wattr_on(win,attrs,opts) result (ires) bind(C, name='wattr_on')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function wattr_on
!         int wattr_set(WINDOW *win, attr_t attrs, short color_pair, void *opts);
    function wattr_set(win,attrs,color_pair,opts) result (ires) bind(C, name='wattr_set')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function wattr_set

!         int chgat(int n, attr_t attr, short color, const void *opts);
    function chgat(n,attr,color,opts) result (ires) bind(C, name='chgat')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function chgat
!         int mvchgat(int y, int x, int n, attr_t attr, short color, const void *opts);
    function mvchgat(y,x,n,attr,color,opts) result (ires) bind(C, name='mvchgat')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function mvchgat
!         int mvwchgat(WINDOW *win, int y, int x, int n, attr_t attr, short color, const void *opts);
    function mvwchgat(win,y,x,n,attr,color,opts) result (ires) bind(C, name='mvwchgat')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function mvwchgat
!         int wchgat(WINDOW *win, int n, attr_t attr, short color, const void *opts);
    function wchgat(win,n,attr,color,opts) result (ires) bind(C, name='wchgat')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function wchgat

!         chtype getattrs(WINDOW *win);
    function getattrs(win) result (ires) bind(C, name='getattrs')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getattrs
!         int     getbegx(WINDOW *);
    function getbegx(win) result (ires) bind(C, name='getbegx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getbegx
!         int     getbegy(WINDOW *);
    function getbegy(win) result (ires) bind(C, name='getbegy')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getbegy
!         int     getmaxx(WINDOW *);
    function getmaxx(win) result (ires) bind(C, name='getmaxx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getmaxx
!         int     getmaxy(WINDOW *);
    function getmaxy(win) result (ires) bind(C, name='getmaxy')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getmaxy
!         int     getparx(WINDOW *);
    function getparx(win) result (ires) bind(C, name='getparx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getparx
!         int     getpary(WINDOW *);
    function getpary(win) result (ires) bind(C, name='getpary')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getpary
!         int     getcurx(WINDOW *);
    function getcurx(win) result (ires) bind(C, name='getcurx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getcurx
!         int     getcury(WINDOW *);
    function getcury(win) result (ires) bind(C, name='getcury')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getcury
!         int beep(void);
    function beep() result (ires) bind(C, name='beep')
      use iso_c_binding
      integer (C_INT) :: ires
    end function beep
!         int flash(void);
    function flash() result (ires) bind(C, name='flash')
      use iso_c_binding
      integer (C_INT) :: ires
    end function flash
!         int bkgd(chtype ch);
    function bkgd(ch) result (ires) bind(C, name='bkgd')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function bkgd
!         void bkgdset(chtype ch);
    subroutine bkgdset(ch) bind(C, name='bkgdset')
      use iso_c_binding
      integer (C_LONG), value :: ch
    end subroutine bkgdset
!         chtype getbkgd(WINDOW *win);
    function getbkgd(win) result (ires) bind(C, name='getbkgd')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getbkgd
!         int wbkgd(WINDOW *win, chtype ch);
    function wbkgd(win,ch) result (ires) bind(C, name='wbkgd')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function wbkgd
!         void wbkgdset(WINDOW *win, chtype ch);
    subroutine wbkgdset(win,ch) bind(C, name='wbkgdset')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end subroutine wbkgdset
!         int box(WINDOW *win, chtype verch, chtype horch);
    function box(win,verch,horch) result (ires) bind(C, name='box')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: verch
      integer (C_LONG), value :: horch
    end function box
!         int hline(chtype ch, int n);
    function hline(ch,n) result (ires) bind(C, name='hline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function hline
!         int vline(chtype ch, int n);
    function vline(ch,n) result (ires) bind(C, name='vline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function vline
!         int whline(WINDOW *win, chtype ch, int n);
    function whline(win,ch,n) result (ires) bind(C, name='whline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function whline
!         int wvline(WINDOW *win, chtype ch, int n);
    function wvline(win,ch,n) result (ires) bind(C, name='wvline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function wvline
!         int mvhline(int y, int x, chtype ch, int n);
    function mvhline(y,x,ch,n) result (ires) bind(C, name='mvhline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvhline
!         int mvvline(int y, int x, chtype ch, int n);
    function mvvline(y,x,ch,n) result (ires) bind(C, name='mvvline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvvline
!         int mvwhline(WINDOW *win, int y, int x, chtype ch, int n);
    function mvwhline(win,y,x,ch,n) result (ires) bind(C, name='mvwhline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvwhline
!         int mvwvline(WINDOW *win, int y, int x, chtype ch, int n);
    function mvwvline(win,y,x,ch,n) result (ires) bind(C, name='mvwvline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvwvline
!         int clear(void);
    function clear() result (ires) bind(C, name='clear')
      use iso_c_binding
      integer (C_INT) :: ires
    end function clear
!         int wclear(WINDOW *win);
    function wclear(win) result (ires) bind(C, name='wclear')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wclear
!         int erase(void);
    function erase() result (ires) bind(C, name='erase')
      use iso_c_binding
      integer (C_INT) :: ires
    end function erase
!         int werase(WINDOW *win);
    function werase(win) result (ires) bind(C, name='werase')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function werase
!         int clrtobot(void);
    function clrtobot() result (ires) bind(C, name='clrtobot')
      use iso_c_binding
      integer (C_INT) :: ires
    end function clrtobot
!         int wclrtobot(WINDOW *win);
    function wclrtobot(win) result (ires) bind(C, name='wclrtobot')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wclrtobot
!         int clrtoeol(void);
    function clrtoeol() result (ires) bind(C, name='clrtoeol')
      use iso_c_binding
      integer (C_INT) :: ires
    end function clrtoeol
!         int wclrtoeol(WINDOW *win);
    function wclrtoeol(win) result (ires) bind(C, name='wclrtoeol')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wclrtoeol
!         int start_color(void);
    function start_color() result (ires) bind(C, name='start_color')
      use iso_c_binding
      integer (C_INT) :: ires
    end function start_color
!         int init_pair(short pair, short fg, short bg);
    function init_pair(pair,fg,bg) result (ires) bind(C, name='init_pair')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: pair
      integer (C_SHORT), value :: fg
      integer (C_SHORT), value :: bg
    end function init_pair
!         int init_color(short color, short red, short green, short blue);
    function init_color(color,red,green,blue) result (ires) bind(C, name='init_color')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color
      integer (C_SHORT), value :: red
      integer (C_SHORT), value :: green
      integer (C_SHORT), value :: blue
    end function init_color
!         bool has_colors(void);
    function has_colors() result (ires) bind(C, name='has_colors')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function has_colors
!         bool can_change_color(void);
    function can_change_color() result (ires) bind(C, name='can_change_color')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function can_change_color
!         int color_content(short color, short *red, short *green, short *blue);
    function color_content(color,red,green,blue) result (ires) bind(C, name='color_content')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color
      type (C_PTR), value :: red
      type (C_PTR), value :: green
      type (C_PTR), value :: blue
    end function color_content
!         int pair_content(short pair, short *fg, short *bg);
    function pair_content(pair,fg,bg) result (ires) bind(C, name='pair_content')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: pair
      type (C_PTR), value :: fg
      type (C_PTR), value :: bg
    end function pair_content

!         int assume_default_colors(int f, int b);
    function assume_default_colors(f,b) result (ires) bind(C, name='assume_default_colors')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: f
      integer (C_INT), value :: b
    end function assume_default_colors
!         int use_default_colors(void);
    function use_default_colors() result (ires) bind(C, name='use_default_colors')
      use iso_c_binding
      integer (C_INT) :: ires
    end function use_default_colors

!         int PDC_set_line_color(short color);
    function PDC_set_line_color(color) result (ires) bind(C, name='PDC_set_line_color')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color
    end function PDC_set_line_color
!         int delch(void);
    function delch() result (ires) bind(C, name='delch')
      use iso_c_binding
      integer (C_INT) :: ires
    end function delch
!         int wdelch(WINDOW *win);
    function wdelch(win) result (ires) bind(C, name='wdelch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wdelch
!         int mvdelch(int y, int x);
    function mvdelch(y,x) result (ires) bind(C, name='mvdelch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvdelch
!         int mvwdelch(WINDOW *win, int y, int x);
    function mvwdelch(win,y,x) result (ires) bind(C, name='mvwdelch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwdelch
!         int deleteln(void);
    function deleteln() result (ires) bind(C, name='deleteln')
      use iso_c_binding
      integer (C_INT) :: ires
    end function deleteln
!         int wdeleteln(WINDOW *win);
    function wdeleteln(win) result (ires) bind(C, name='wdeleteln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wdeleteln
!         int insdelln(int n);
    function insdelln(n) result (ires) bind(C, name='insdelln')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: n
    end function insdelln
!         int winsdelln(WINDOW *win, int n);
    function winsdelln(win,n) result (ires) bind(C, name='winsdelln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: n
    end function winsdelln
!         int insertln(void);
    function insertln() result (ires) bind(C, name='insertln')
      use iso_c_binding
      integer (C_INT) :: ires
    end function insertln
!         int winsertln(WINDOW *win);
    function winsertln(win) result (ires) bind(C, name='winsertln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function winsertln

!         int mvdeleteln(int y, int x);
    function mvdeleteln(y,x) result (ires) bind(C, name='mvdeleteln')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvdeleteln
!         int mvwdeleteln(WINDOW *win, int y, int x);
    function mvwdeleteln(win,y,x) result (ires) bind(C, name='mvwdeleteln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwdeleteln
!         int mvinsertln(int y, int x);
    function mvinsertln(y,x) result (ires) bind(C, name='mvinsertln')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvinsertln
!         int mvwinsertln(WINDOW *win, int y, int x);
    function mvwinsertln(win,y,x) result (ires) bind(C, name='mvwinsertln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwinsertln
!         int wgetch(WINDOW *win);
    function wgetch(win) result (ires) bind(C, name='wgetch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wgetch
!         int mvgetch(int y, int x);
    function mvgetch(y,x) result (ires) bind(C, name='mvgetch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvgetch
!         int mvwgetch(WINDOW *win, int y, int x);
    function mvwgetch(win,y,x) result (ires) bind(C, name='mvwgetch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwgetch
!         int ungetch(int ch);
    function ungetch(ch) result (ires) bind(C, name='ungetch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ch
    end function ungetch
!         int flushinp(void);
    function flushinp() result (ires) bind(C, name='flushinp')
      use iso_c_binding
      integer (C_INT) :: ires
    end function flushinp

!         int getstr(char *str);
    function getstr(str) result (ires) bind(C, name='getstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
    end function getstr
!         int wgetstr(WINDOW *win, char *str);
    function wgetstr(win,str) result (ires) bind(C, name='wgetstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
    end function wgetstr
!         int mvgetstr(int y, int x, char *str);
    function mvgetstr(y,x,str) result (ires) bind(C, name='mvgetstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvgetstr
!         int mvwgetstr(WINDOW *win, int y, int x, char *str);
    function mvwgetstr(win,y,x,str) result (ires) bind(C, name='mvwgetstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvwgetstr
!         int getnstr(char *str, int n);
    function getnstr(str,n) result (ires) bind(C, name='getnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function getnstr
!         int wgetnstr(WINDOW *win, char *str, int n);
    function wgetnstr(win,str,n) result (ires) bind(C, name='wgetnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function wgetnstr
!         int mvgetnstr(int y, int x, char *str, int n);
    function mvgetnstr(y,x,str,n) result (ires) bind(C, name='mvgetnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvgetnstr
!         int mvwgetnstr(WINDOW *win, int y, int x, char *str, int n);
    function mvwgetnstr(win,y,x,str,n) result (ires) bind(C, name='mvwgetnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvwgetnstr
!
! needs to be defined (is macro)
!
!         void getyx(WINDOW *win, int y, int x);
    subroutine getyx(win,y,x) bind(C, name='getyx')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end subroutine getyx

!         chtype inch(void);
    function inch() result (ires) bind(C, name='inch')
      use iso_c_binding
      integer (C_LONG) :: ires
    end function inch
!         chtype winch(WINDOW *win);
    function winch(win) result (ires) bind(C, name='winch')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function winch
!         chtype mvinch(int y, int x);
    function mvinch(y,x) result (ires) bind(C, name='mvinch')
      use iso_c_binding
      integer (C_LONG) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvinch
!         chtype mvwinch(WINDOW *win, int y, int x);
    function mvwinch(win,y,x) result (ires) bind(C, name='mvwinch')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwinch

!         int inchstr(chtype *ch);
    function inchstr(ch) result (ires) bind(C, name='inchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: ch
    end function inchstr
!         int inchnstr(chtype *ch, int n);
    function inchnstr(ch,n) result (ires) bind(C, name='inchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function inchnstr
!         int winchstr(WINDOW *win, chtype *ch);
    function winchstr(win,ch) result (ires) bind(C, name='winchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), value :: ch
    end function winchstr
!         int winchnstr(WINDOW *win, chtype *ch, int n);
    function winchnstr(win,ch,n) result (ires) bind(C, name='winchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function winchnstr
!         int mvinchstr(int y, int x, chtype *ch);
    function mvinchstr(y,x,ch) result (ires) bind(C, name='mvinchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
    end function mvinchstr
!         int mvinchnstr(int y, int x, chtype *ch, int n);
    function mvinchnstr(y,x,ch,n) result (ires) bind(C, name='mvinchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function mvinchnstr
!         int mvwinchstr(WINDOW *, int y, int x, chtype *ch);
    function mvwinchstr(win,y,x,ch) result (ires) bind(C, name='mvwinchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
    end function mvwinchstr
!         int mvwinchnstr(WINDOW *, int y, int x, chtype *ch, int n);
    function mvwinchnstr(win,y,x,ch,n) result (ires) bind(C, name='mvwinchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function mvwinchnstr
!         WINDOW *initscr(void);
    function initscr() result (iwin) bind(C, name='initscr')
      use iso_c_binding
      type (C_PTR) :: iwin
    end function initscr
!         WINDOW *Xinitscr(int argc, char *argv[]);
!   function Xinitscr(argc,argv) result (iwin) bind(C, name='Xinitscr')
!     use iso_c_binding
!     type (C_PTR) :: iwin
!     integer (C_INT), value :: argc
!     character (C_CHAR), dimension(*) :: argv
!   end function Xinitscr
!         int endwin(void);
    function endwin() result (ires) bind(C, name='endwin')
      use iso_c_binding
      integer (C_INT) :: ires
    end function endwin
!         bool isendwin(void);
    function isendwin() result (ires) bind(C, name='isendwin')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function isendwin
!         int resize_term(int nlines, int ncols);
    function resize_term(nlines,ncols) result (ires) bind(C, name='resize_term')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
    end function resize_term
!         bool is_termresized(void);
    function is_termresized() result (ires) bind(C, name='is_termresized')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function is_termresized
!         int cbreak(void);
    function cbreak() result (ires) bind(C, name='cbreak')
      use iso_c_binding
      integer (C_INT) :: ires
    end function cbreak
!         int nocbreak(void);
    function nocbreak() result (ires) bind(C, name='nocbreak')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nocbreak
!         int echo(void);
    function echo() result (ires) bind(C, name='echo')
      use iso_c_binding
      integer (C_INT) :: ires
    end function echo
!         int noecho(void);
    function noecho() result (ires) bind(C, name='noecho')
      use iso_c_binding
      integer (C_INT) :: ires
    end function noecho
!         int halfdelay(int tenths);
    function halfdelay(tenths) result (ires) bind(C, name='halfdelay')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: tenths
    end function halfdelay
!         int intrflush(WINDOW *win, bool bf);
    function intrflush(win, bf) result (ires) bind(C, name='intrflush')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function intrflush
!         int keypad(WINDOW *win, bool bf);
    function keypad(win, bf) result (ires) bind(C, name='keypad')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function keypad
!         int meta(WINDOW *win, bool bf);
    function meta(win, bf) result (ires) bind(C, name='meta')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function meta
!         int nl(void);
    function nl() result (ires) bind(C, name='nl')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nl
!         int nonl(void);
    function nonl() result (ires) bind(C, name='nonl')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nonl
!         int nodelay(WINDOW *win, bool bf);
    function nodelay(win, bf) result (ires) bind(C, name='nodelay')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function nodelay
!         int notimeout(WINDOW *win, bool bf);
    function notimeout(win, bf) result (ires) bind(C, name='notimeout')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function notimeout
!         int raw(void);
    function raw() result (ires) bind(C, name='raw')
      use iso_c_binding
      integer (C_INT) :: ires
    end function raw
!         int noraw(void);
    function noraw() result (ires) bind(C, name='noraw')
      use iso_c_binding
      integer (C_INT) :: ires
    end function noraw
!         void noqiflush(void);
    subroutine noqiflush() bind(C, name='noqiflush')
      use iso_c_binding
    end subroutine noqiflush
!         void qiflush(void);
    subroutine qiflush() bind(C, name='qiflush')
      use iso_c_binding
    end subroutine qiflush
!         void timeout(int delay);
    subroutine timeout(delay) bind(C, name='timeout')
      use iso_c_binding
      integer (C_INT), value :: delay
    end subroutine timeout
!         void wtimeout(WINDOW *win, int delay);
    subroutine wtimeout(win,delay) bind(C, name='wtimeout')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_INT), value :: delay
    end subroutine wtimeout
!         int typeahead(int fildes);
    function typeahead(fildes) result (ires) bind(C, name='typeahead')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: fildes
    end function typeahead

!         int crmode(void);
    function crmode() result (ires) bind(C, name='crmode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function crmode
!         int nocrmode(void);
    function nocrmode() result (ires) bind(C, name='nocrmode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nocrmode
!         int insch(chtype ch);
    function insch(ch) result (ires) bind(C, name='insch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function insch
!         int winsch(WINDOW *win, chtype ch);
    function winsch(win,ch) result (ires) bind(C, name='winsch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function winsch
!         int mvinsch(int y, int x, chtype ch);
    function mvinsch(y,x,ch) result (ires) bind(C, name='mvinsch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvinsch
!         int mvwinsch(WINDOW *win, int y, int x, chtype ch);
    function mvwinsch(win,y,x,ch) result (ires) bind(C, name='mvwinsch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvwinsch

!         int insrawch(chtype ch);
    function insrawch(ch) result (ires) bind(C, name='insrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function insrawch
!         int winsrawch(WINDOW *win, chtype ch);
    function winsrawch(win,ch) result (ires) bind(C, name='winsrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function winsrawch
!         int mvinsrawch(int y, int x, chtype ch);
    function mvinsrawch(y,x,ch) result (ires) bind(C, name='mvinsrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvinsrawch
!         int mvwinsrawch(WINDOW *win, int y, int x, chtype ch);
    function mvwinsrawch(win,y,x,ch) result (ires) bind(C, name='mvwinsrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvwinsrawch
!         int insstr(const char *str);
    function insstr(str) result (ires) bind(C, name='insstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
    end function insstr
!         int insnstr(const char *str, int n);
    function insnstr(str,n) result (ires) bind(C, name='insnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function insnstr
!         int winsstr(WINDOW *win, const char *str);
    function winsstr(win,str) result (ires) bind(C, name='winsstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
    end function winsstr
!         int winsnstr(WINDOW *win, const char *str, int n);
    function winsnstr(win,str,n) result (ires) bind(C, name='winsnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function winsnstr
!         int mvinsstr(int y, int x, const char *str);
    function mvinsstr(y,x,str) result (ires) bind(C, name='mvinsstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvinsstr
!         int mvinsnstr(int y, int x, const char *str, int n);
    function mvinsnstr(y,x,str,n) result (ires) bind(C, name='mvinsnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvinsnstr
!         int mvwinsstr(WINDOW *win, int y, int x, const char *str);
    function mvwinsstr(win,y,x,str) result (ires) bind(C, name='mvwinsstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvwinsstr
!         int mvwinsnstr(WINDOW *win, int y, int x, const char *str, int n);
    function mvwinsnstr(win,y,x,str,n) result (ires) bind(C, name='mvwinsnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvwinsnstr

!         int instr(char *str);
    function instr(str) result (ires) bind(C, name='instr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
    end function instr
!         int innstr(char *str, int n);
    function innstr(str,n) result (ires) bind(C, name='innstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function innstr
!         int winstr(WINDOW *win, char *str);
    function winstr(win,str) result (ires) bind(C, name='winstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
    end function winstr
!         int winnstr(WINDOW *win, char *str, int n);
    function winnstr(win,str,n) result (ires) bind(C, name='winnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function winnstr
!         int mvinstr(int y, int x, char *str);
    function mvinstr(y,x,str) result (ires) bind(C, name='mvinstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvinstr
!         int mvinnstr(int y, int x, char *str, int n);
    function mvinnstr(y,x,str,n) result (ires) bind(C, name='mvinnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvinnstr
!         int mvwinstr(WINDOW *win, int y, int x, char *str);
    function mvwinstr(win,y,x,str) result (ires) bind(C, name='mvwinstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvwinstr
!         int mvwinnstr(WINDOW *win, int y, int x, char *str, int n);
    function mvwinnstr(win,y,x,str,n) result (ires) bind(C, name='mvwinnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvwinnstr
!         int def_prog_mode(void);
    function def_prog_mode() result (ires) bind(C, name='def_prog_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function def_prog_mode
!         int def_shell_mode(void);
    function def_shell_mode() result (ires) bind(C, name='def_shell_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function def_shell_mode
!         int reset_prog_mode(void);
    function reset_prog_mode() result (ires) bind(C, name='reset_prog_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function reset_prog_mode
!         int reset_shell_mode(void);
    function reset_shell_mode() result (ires) bind(C, name='reset_shell_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function reset_shell_mode
!         int resetty(void);
    function resetty() result (ires) bind(C, name='resetty')
      use iso_c_binding
      integer (C_INT) :: ires
    end function resetty
!         int savetty(void);
    function savetty() result (ires) bind(C, name='savetty')
      use iso_c_binding
      integer (C_INT) :: ires
    end function savetty
!         int ripoffline(int line, int (*init)(WINDOW *, int));
!         int copywin(const WINDOW *, WINDOW *, int, int, int, int, int, int, int);
    function copywin(srcwin, dstwin, sminrow, smincol, dminrow, dmincol,  &
                     dmaxrow, dmaxcol, overlay) result (ires) bind(C, name='copywin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value, intent(in) :: srcwin
      type (C_PTR), value :: dstwin
      integer (C_INT), value :: sminrow, smincol
      integer (C_INT), value :: dminrow, dmincol
      integer (C_INT), value :: dmaxrow, dmaxcol
      integer (C_INT), value :: overlay
    end function copywin
!         int overlay(const WINDOW *srcwin, WINDOW *dstwin);
    function overlay(srcwin, dstwin) result (ires) bind(C, name='overlay')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value, intent(in) :: srcwin
      type (C_PTR), value :: dstwin
    end function overlay
!         int overwrite(const WINDOW *srcwin, WINDOW *dstwin);
    function overwrite(srcwin, dstwin) result (ires) bind(C, name='overwrite')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value, intent(in) :: srcwin
      type (C_PTR), value :: dstwin
    end function overwrite
!         int curs_set(int visibility);
    function curs_set(visibility) result (ires) bind(C, name='curs_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: visibility
    end function curs_set
!         int napms(int ms);
    function napms(ms) result (ires) bind(C, name='napms')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ms
    end function napms

!         int draino(int ms);
    function draino(ms) result (ires) bind(C, name='draino')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ms
    end function draino
!         int resetterm(void);
    function resetterm() result (ires) bind(C, name='resetterm')
      use iso_c_binding
      integer (C_INT) :: ires
    end function resetterm
!         int fixterm(void);
    function fixterm() result (ires) bind(C, name='fixterm')
      use iso_c_binding
      integer (C_INT) :: ires
    end function fixterm
!         int saveterm(void);
    function saveterm() result (ires) bind(C, name='saveterm')
      use iso_c_binding
      integer (C_INT) :: ires
    end function saveterm
!         char *keyname(int key);
    function keyname(key) result (ich) bind(C, name='keyname')
      use iso_c_binding
      character (C_CHAR) :: ich
      integer (C_INT), value :: key
    end function keyname
!         bool has_key(int key);
    function has_key(key) result (ires) bind(C, name='has_key')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      integer (C_INT), value :: key
    end function has_key
!         int mouse_set(unsigned long mbe);
    function mouse_set(mbe) result (ires) bind(C, name='mouse_set')
      use iso_c_binding
      integer (C_LONG) :: mbe
      integer (C_INT) :: ires
    end function mouse_set
!         int mouse_on(unsigned long mbe);
    function mouse_on(mbe) result (ires) bind(C, name='mouse_on')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG) :: mbe
    end function mouse_on
!         int mouse_off(unsigned long mbe);
    function mouse_off(mbe) result (ires) bind(C, name='mouse_off')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG) :: mbe
    end function mouse_off
!         int request_mouse_pos(void);
    function request_mouse_pos() result (ires) bind(C, name='request_mouse_pos')
      use iso_c_binding
      integer (C_INT) :: ires
    end function request_mouse_pos
!         int map_button(unsigned long button);
    function map_button(button) result (ires) bind(C, name='map_button')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG) :: button
    end function map_button
!         void wmouse_position(WINDOW *win, int *y, int *x);
    subroutine wmouse_position(win,y,x) bind(C, name='wmouse_position')
      use iso_c_binding
      type (C_PTR), value :: win
      type (C_PTR), value :: y
      type (C_PTR), value :: x
    end subroutine wmouse_position
!         unsigned long getmouse(void);
!         bool wenclose(const WINDOW *win, int y, int x);
    function wenclose(win,y,x) result (ires) bind(C, name='wenclose')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      type (C_PTR), intent(in) :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function wenclose
!         int move(int y, int x);
    function move(y,x) result (ires) bind(C, name='move')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function move
!         int wmove(WINDOW *win, int y, int x);
    function wmove(win,y,x) result (ires) bind(C, name='wmove')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function wmove

!         int clearok(WINDOW *win, bool bf);
    function clearok(win, bf) result (ires) bind(C, name='clearok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function clearok
!         int idlok(WINDOW *win, bool bf);
    function idlok(win, bf) result (ires) bind(C, name='idlok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function idlok
!         void idcok(WINDOW *win, bool bf);
    subroutine idcok(win, bf) bind(C, name='idcok')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end subroutine idcok
!         void immedok(WINDOW *win, bool bf);
    subroutine immedok(win, bf) bind(C, name='immedok')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end subroutine immedok
!         int leaveok(WINDOW *win, bool bf);
    function leaveok(win, bf) result (ires) bind(C, name='leaveok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function leaveok
!         int setscrreg(int top, int bot);
    function setscrreg(top,bot) result (ires) bind(C, name='setscrreg')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: top
      integer (C_INT), value :: bot
    end function setscrreg
!         int wsetscrreg(WINDOW *win, int top, int bot);
    function wsetscrreg(win,top,bot) result (ires) bind(C, name='wsetscrreg')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: top
      integer (C_INT), value :: bot
    end function wsetscrreg
!         int scrollok(WINDOW *win, bool bf);
    function scrollok(win, bf) result (ires) bind(C, name='scrollok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function scrollok

!         int raw_output(bool bf);
    function raw_output(bf) result (ires) bind(C, name='raw_output')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SIGNED_CHAR), value :: bf
    end function raw_output
!         int refresh(void);
    function refresh() result (ires) bind(C, name='refresh')
      use iso_c_binding
      integer (C_INT) :: ires
    end function refresh
!         int wrefresh(WINDOW *win);
    function wrefresh(win) result (ires) bind(C, name='wrefresh')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wrefresh
!         int wnoutrefresh(WINDOW *win);
    function wnoutrefresh(win) result (ires) bind(C, name='wnoutrefresh')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wnoutrefresh
!         int doupdate(void);
    function doupdate() result (ires) bind(C, name='doupdate')
      use iso_c_binding
      integer (C_INT) :: ires
    end function doupdate
!         int redrawwin(WINDOW *win);
    function redrawwin(win) result (ires) bind(C, name='redrawwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function redrawwin
!         int wredrawln(WINDOW *win, int beg_line, int num_lines);
    function wredrawln(win,beg_line,num_lines) result (ires) bind(C, name='wredrawln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: beg_line
      integer (C_INT), value :: num_lines
    end function wredrawln

!         int scanw(const char *fmt, ...);
    function scanw(formt) result (ires) bind(C, name='scanw')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function scanw
!         int wscanw(WINDOW *win, const char *fmt, ...);
    function wscanw(win,formt) result (ires) bind(C, name='wscanw')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function wscanw
!         int mvscanw(int y, int x, const char *fmt, ...);
    function mvscanw(y,x,formt) result (ires) bind(C, name='mvscanw')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function mvscanw
!         int mvwscanw(WINDOW *win, int y, int x, const char *fmt, ...);
    function mvwscanw(win,y,x,formt) result (ires) bind(C, name='mvwscanw')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function mvwscanw
!         int vwscanw(WINDOW *win, const char *fmt, va_list varglist);
!         int vw_scanw(WINDOW *win, const char *fmt, va_list varglist);
!         int scroll(WINDOW *win);
    function scroll(win) result (ires) bind(C, name='scroll')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function scroll
!         int scrl(int n);
    function scrl(n) result (ires) bind(C, name='scrl')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: n
    end function scrl
!         int wscrl(WINDOW *win, int n);
    function wscrl(win,n) result (ires) bind(C, name='wscrl')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: n
    end function wscrl

!         void filter(void);
    subroutine filter() bind(C, name='filter')
      use iso_c_binding
    end subroutine filter
!         void use_env(bool x);
    subroutine use_env(x) bind(C, name='use_env')
      use iso_c_binding
      integer (C_SIGNED_CHAR), value :: x
    end subroutine use_env
!         int delay_output(int ms);
    function delay_output(ms) result (ires) bind(C, name='delay_output')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ms
    end function delay_output
!         int slk_init(int fmt);
    function slk_init(fmt) result (ires) bind(C, name='slk_init')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: fmt
    end function slk_init
!         int slk_set(int labnum, const char *label, int justify);
    function slk_set(labnum,label,justify) result (ires) bind(C, name='slk_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: labnum
      character (C_CHAR), dimension(*), intent(in) :: label
      integer (C_INT), value :: justify
    end function slk_set
!         int slk_refresh(void);
    function slk_refresh() result (ires) bind(C, name='slk_refresh')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_refresh
!         int slk_noutrefresh(void);
    function slk_noutrefresh() result (ires) bind(C, name='slk_noutrefresh')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_noutrefresh
!         char *slk_label(int labnum);
!         int slk_clear(void);
    function slk_clear() result (ires) bind(C, name='slk_clear')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_clear
!         int slk_restore(void);
    function slk_restore() result (ires) bind(C, name='slk_restore')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_restore
!         int slk_touch(void);
    function slk_touch() result (ires) bind(C, name='slk_touch')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_touch
!         int slk_attron(const chtype attrs);
    function slk_attron(attrs) result (ires) bind(C, name='slk_attron')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
    end function slk_attron
!         int slk_attr_on(const attr_t attrs, void *opts);
    function slk_attr_on(attrs,opts) result (ires) bind(C, name='slk_attr_on')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
      type (C_PTR), value :: opts
    end function slk_attr_on
!         int slk_attrset(const chtype attrs);
    function slk_attrset(attrs) result (ires) bind(C, name='slk_attrset')
      use iso_c_binding
      integer (C_INT) :: ires
      integ
module curses
  use iso_c_binding
  implicit none

  integer (C_SIGNED_CHAR), parameter :: curses_true = 1
  integer (C_SIGNED_CHAR), parameter :: curses_false = 0
  integer, parameter :: curses_err = -1
  integer, parameter :: curses_ok  =  0
! ---------------------------------------------------------------------
!
!  PDCurses Mouse Interface -- SYSVR4, with extensions
!
  type mouse_data
    integer(C_INT) :: x        ! absolute column, 0 based, measured in characters
    integer(C_INT) :: y        ! absolute row, 0 based, measured in characters
    integer(C_SHORT), dimension(3) :: button ! state of each button
    integer(C_INT) :: changes  ! flags indicating what has changed with the mouse
  end type mouse_data

  type (mouse_data) :: mouse_status

!------------------------------------------------------------------------
!
!   PDCurses Structure Definitions
!
! definition of a window
! type window_data
!   integer (C_INT) :: cury, curx    ! current pseudo-cursor
!   integer (C_INT) :: maxy, maxx    ! max window coordinates
!   integer (C_INT) :: begy, begx    ! origin on screen
!   integer (C_INT) :: flags         ! window properties
!   integer (C_LONG) ::  attrs       ! standard attributes and colors
!   integer (C_LONG) ::  bkgd        ! background, normally blank
!   integer (C_SIGNED_CHAR) :: clear      ! causes clear at next refresh
!   integer (C_SIGNED_CHAR) :: leaveit    ! leaves cursor where it is
!   integer (C_SIGNED_CHAR) :: scroll     ! allows window scrolling
!   integer (C_SIGNED_CHAR) :: nodelay    ! input character wait flag
!   integer (C_SIGNED_CHAR) :: immed      ! immediate update flag
!   integer (C_SIGNED_CHAR) :: sync       ! synchronise window ancestors
!   integer (C_SIGNED_CHAR) :: use_keypad ! flags keypad key mode active
!   type (C_PTR) :: y                ! pointer to line pointer array
!   type (C_PTR) :: firstch          ! first changed character in line
!   type (C_PTR) :: lastch           ! last changed character in line
!   integer (C_INT) :: tmarg         ! top of scrolling region
!   integer (C_INT) :: bmarg         ! bottom of scrolling region
!   integer (C_INT) :: delayms       ! milliseconds of delay for getch()
!   integer (C_INT) :: parx, pary    ! coords relative to parent (0,0)
!   type (C_PTR) :: parent           ! subwin pointer to parent win
! end type window_data
!
  integer, parameter :: curses_button_released         =  0
  integer, parameter :: curses_button_pressed          =  1
  integer, parameter :: curses_button_clicked          =  2
  integer, parameter :: curses_button_double_clicked   =  3
  integer, parameter :: curses_button_triple_clicked   =  4
  integer, parameter :: curses_button_moved            =  5
  integer, parameter :: curses_wheel_scrolled          =  6
  integer, parameter :: curses_button_action_mask      =  7
  integer, parameter :: curses_pdc_button_shift        =  8
  integer, parameter :: curses_pdc_button_control      = 16
  integer, parameter :: curses_pdc_button_alt          = 32
  integer, parameter :: curses_button_modifier_mask    = 56

  integer, parameter :: curses_color_black     = 0
  integer, parameter :: curses_color_red       = 1
  integer, parameter :: curses_color_green     = 2
  integer, parameter :: curses_color_yellow    = 3
  integer, parameter :: curses_color_blue      = 4
  integer, parameter :: curses_color_magenta   = 5
  integer, parameter :: curses_color_cyan      = 6
  integer, parameter :: curses_color_white     = 7

  integer, parameter :: curses_key_code_yes  =  256 !  If get_wch() gives a key code

  integer, parameter :: curses_key_break     =  257 ! Not on PC KBD
  integer, parameter :: curses_key_down      =  258 ! Down arrow key
  integer, parameter :: curses_key_up        =  259 ! Up arrow key
  integer, parameter :: curses_key_left      =  260 ! Left arrow key
  integer, parameter :: curses_key_right     =  261 ! Right arrow key
  integer, parameter :: curses_key_home      =  262 ! home key
  integer, parameter :: curses_key_backspace =  263 ! not on pc
  integer, parameter :: curses_key_f0        =  264 ! function keys; 64 reserved

  integer, parameter :: curses_key_dl        =  328 ! delete line
  integer, parameter :: curses_key_il        =  329 ! insert line
  integer, parameter :: curses_key_dc        =  330 ! delete character
  integer, parameter :: curses_key_ic        =  331 ! insert char or enter ins mode
  integer, parameter :: curses_key_eic       =  332 ! exit insert char mode
  integer, parameter :: curses_key_clear     =  333 ! clear screen
  integer, parameter :: curses_key_eos       =  334 ! clear to end of screen
  integer, parameter :: curses_key_eol       =  335 ! clear to end of line
  integer, parameter :: curses_key_sf        =  336 ! scroll 1 line forward
  integer, parameter :: curses_key_sr        =  337 ! scroll 1 line back (reverse)
  integer, parameter :: curses_key_npage     =  338 ! next page
  integer, parameter :: curses_key_ppage     =  339 ! previous page
  integer, parameter :: curses_key_stab      =  340 ! set tab
  integer, parameter :: curses_key_ctab      =  341 ! clear tab
  integer, parameter :: curses_key_catab     =  342 ! clear all tabs
  integer, parameter :: curses_key_enter     =  343 ! enter or send (unreliable)
  integer, parameter :: curses_key_sreset    =  344 ! soft/reset (partial/unreliable)
  integer, parameter :: curses_key_reset     =  345 ! reset/hard reset (unreliable)
  integer, parameter :: curses_key_print     =  346 ! print/copy
  integer, parameter :: curses_key_ll        =  347 ! home down/bottom (lower left)
  integer, parameter :: curses_key_abort     =  348 ! abort/terminate key (any)
  integer, parameter :: curses_key_shelp     =  349 ! short help
  integer, parameter :: curses_key_lhelp     =  350 ! long help
  integer, parameter :: curses_key_btab      =  351 ! Back tab key
  integer, parameter :: curses_key_beg       =  352 ! beg(inning) key
  integer, parameter :: curses_key_cancel    =  353 ! cancel key
  integer, parameter :: curses_key_close     =  354 ! close key
  integer, parameter :: curses_key_command   =  355 ! cmd (command) key
  integer, parameter :: curses_key_copy      =  356 ! copy key
  integer, parameter :: curses_key_create    =  357 ! create key
  integer, parameter :: curses_key_end       =  358 ! end key
  integer, parameter :: curses_key_exit      =  359 ! exit key
  integer, parameter :: curses_key_find      =  360 ! find key
  integer, parameter :: curses_key_help      =  361 ! help key
  integer, parameter :: curses_key_mark      =  362 ! mark key
  integer, parameter :: curses_key_message   =  363 ! message key
  integer, parameter :: curses_key_move      =  364 ! move key
  integer, parameter :: curses_key_next      =  365 ! next object key
  integer, parameter :: curses_key_open      =  366 ! open key
  integer, parameter :: curses_key_options   =  367 ! options key
  integer, parameter :: curses_key_previous  =  368 ! previous object key
  integer, parameter :: curses_key_redo      =  369 ! redo key
  integer, parameter :: curses_key_reference =  370 ! ref(erence) key
  integer, parameter :: curses_key_refresh   =  371 ! refresh key
  integer, parameter :: curses_key_replace   =  372 ! replace key
  integer, parameter :: curses_key_restart   =  373 ! restart key
  integer, parameter :: curses_key_resume    =  374 ! resume key
  integer, parameter :: curses_key_save      =  375 ! save key
  integer, parameter :: curses_key_sbeg      =  376 ! shifted beginning key
  integer, parameter :: curses_key_scancel   =  377 ! shifted cancel key
  integer, parameter :: curses_key_scommand  =  378 ! shifted command key
  integer, parameter :: curses_key_scopy     =  379 ! shifted copy key
  integer, parameter :: curses_key_screate   =  380 ! shifted create key
  integer, parameter :: curses_key_sdc       =  381 ! shifted delete char key
  integer, parameter :: curses_key_sdl       =  382 ! shifted delete line key
  integer, parameter :: curses_key_select    =  383 ! select key
  integer, parameter :: curses_key_send      =  384 ! shifted end key
  integer, parameter :: curses_key_seol      =  385 ! shifted clear line key
  integer, parameter :: curses_key_sexit     =  386 ! shifted exit key
  integer, parameter :: curses_key_sfind     =  387 ! shifted find key
  integer, parameter :: curses_key_shome     =  388 ! shifted home key
  integer, parameter :: curses_key_sic       =  389 ! shifted input key

  integer, parameter :: curses_key_sleft     =  391 ! shifted left arrow key
  integer, parameter :: curses_key_smessage  =  392 ! shifted message key
  integer, parameter :: curses_key_smove     =  393 ! shifted move key
  integer, parameter :: curses_key_snext     =  394 ! shifted next key
  integer, parameter :: curses_key_soptions  =  395 ! shifted options key
  integer, parameter :: curses_key_sprevious =  396 ! shifted prev key
  integer, parameter :: curses_key_sprint    =  397 ! shifted print key
  integer, parameter :: curses_key_sredo     =  398 ! shifted redo key
  integer, parameter :: curses_key_sreplace  =  399 ! shifted replace key
  integer, parameter :: curses_key_sright    =  400 ! shifted right arrow
  integer, parameter :: curses_key_srsume    =  401 ! shifted resume key
  integer, parameter :: curses_key_ssave     =  402 ! shifted save key
  integer, parameter :: curses_key_ssuspend  =  403 ! shifted suspend key
  integer, parameter :: curses_key_sundo     =  404 ! shifted undo key
  integer, parameter :: curses_key_suspend   =  405 ! suspend key
  integer, parameter :: curses_key_undo      =  406 ! undo key
!
!  PDCurses-specific key definitions -- PC only
!
  integer, parameter :: curses_alt_0         =  407
  integer, parameter :: curses_alt_1         =  408
  integer, parameter :: curses_alt_2         =  409
  integer, parameter :: curses_alt_3         =  410
  integer, parameter :: curses_alt_4         =  411
  integer, parameter :: curses_alt_5         =  412
  integer, parameter :: curses_alt_6         =  413
  integer, parameter :: curses_alt_7         =  414
  integer, parameter :: curses_alt_8         =  415
  integer, parameter :: curses_alt_9         =  416
  integer, parameter :: curses_alt_a         =  417
  integer, parameter :: curses_alt_b         =  418
  integer, parameter :: curses_alt_c         =  419
  integer, parameter :: curses_alt_d         =  420
  integer, parameter :: curses_alt_e         =  421
  integer, parameter :: curses_alt_f         =  422
  integer, parameter :: curses_alt_g         =  423
  integer, parameter :: curses_alt_h         =  424
  integer, parameter :: curses_alt_i         =  425
  integer, parameter :: curses_alt_j         =  426
  integer, parameter :: curses_alt_k         =  427
  integer, parameter :: curses_alt_l         =  428
  integer, parameter :: curses_alt_m         =  429
  integer, parameter :: curses_alt_n         =  430
  integer, parameter :: curses_alt_o         =  431
  integer, parameter :: curses_alt_p         =  432
  integer, parameter :: curses_alt_q         =  433
  integer, parameter :: curses_alt_r         =  434
  integer, parameter :: curses_alt_s         =  435
  integer, parameter :: curses_alt_t         =  436
  integer, parameter :: curses_alt_u         =  437
  integer, parameter :: curses_alt_v         =  438
  integer, parameter :: curses_alt_w         =  439
  integer, parameter :: curses_alt_x         =  440
  integer, parameter :: curses_alt_y         =  441
  integer, parameter :: curses_alt_z         =  442
  integer, parameter :: curses_ctl_left      =  443 ! Control-Left-Arrow
  integer, parameter :: curses_ctl_right     =  444
  integer, parameter :: curses_ctl_pgup      =  445
  integer, parameter :: curses_ctl_pgdn      =  446
  integer, parameter :: curses_ctl_home      =  447
  integer, parameter :: curses_ctl_end       =  448

  integer, parameter :: curses_key_a1        =  449 ! upper left on Virtual keypad
  integer, parameter :: curses_key_a2        =  450 ! upper middle on Virt. keypad
  integer, parameter :: curses_key_a3        =  451 ! upper right on Vir. keypad
  integer, parameter :: curses_key_b1        =  452 ! middle left on Virt. keypad
  integer, parameter :: curses_key_b2        =  453 ! center on Virt. keypad
  integer, parameter :: curses_key_b3        =  454 ! middle right on Vir. keypad
  integer, parameter :: curses_key_c1        =  455 ! lower left on Virt. keypad
  integer, parameter :: curses_key_c2        =  456 ! lower middle on Virt. keypad
  integer, parameter :: curses_key_c3        =  457 ! lower right on Vir. keypad

  integer, parameter :: curses_padslash      =  458 ! slash on keypad
  integer, parameter :: curses_padenter      =  459 ! enter on keypad
  integer, parameter :: curses_ctl_padenter  =  460 ! ctl-enter on keypad
  integer, parameter :: curses_alt_padenter  =  461 ! alt-enter on keypad
  integer, parameter :: curses_padstop       =  462 ! stop on keypad
  integer, parameter :: curses_padstar       =  463 ! star on keypad
  integer, parameter :: curses_padminus      =  464 ! minus on keypad
  integer, parameter :: curses_padplus       =  465 ! plus on keypad
  integer, parameter :: curses_ctl_padstop   =  466 ! ctl-stop on keypad
  integer, parameter :: curses_ctl_padcenter =  467 ! ctl-enter on keypad
  integer, parameter :: curses_ctl_padplus   =  468 ! ctl-plus on keypad
  integer, parameter :: curses_ctl_padminus  =  469 ! ctl-minus on keypad
  integer, parameter :: curses_ctl_padslash  =  470 ! ctl-slash on keypad
  integer, parameter :: curses_ctl_padstar   =  471 ! ctl-star on keypad
  integer, parameter :: curses_alt_padplus   =  472 ! alt-plus on keypad
  integer, parameter :: curses_alt_padminus  =  473 ! alt-minus on keypad
  integer, parameter :: curses_alt_padslash  =  474 ! alt-slash on keypad
  integer, parameter :: curses_alt_padstar   =  475 ! alt-star on keypad
  integer, parameter :: curses_alt_padstop   =  476 ! alt-stop on keypad
  integer, parameter :: curses_ctl_ins       =  477 ! ctl-insert
  integer, parameter :: curses_alt_del       =  478 ! alt-delete
  integer, parameter :: curses_alt_ins       =  479 ! alt-insert
  integer, parameter :: curses_ctl_up        =  480 ! ctl-up arrow
  integer, parameter :: curses_ctl_down      =  481 ! ctl-down arrow
  integer, parameter :: curses_ctl_tab       =  482 ! ctl-tab
  integer, parameter :: curses_alt_tab       =  483
  integer, parameter :: curses_alt_minus     =  484
  integer, parameter :: curses_alt_equal     =  485
  integer, parameter :: curses_alt_home      =  486
  integer, parameter :: curses_alt_pgup      =  487
  integer, parameter :: curses_alt_pgdn      =  488
  integer, parameter :: curses_alt_end       =  489
  integer, parameter :: curses_alt_up        =  490 ! alt-up arrow
  integer, parameter :: curses_alt_down      =  491 ! alt-down arrow
  integer, parameter :: curses_alt_right     =  492 ! alt-right arrow
  integer, parameter :: curses_alt_left      =  493 ! alt-left arrow
  integer, parameter :: curses_alt_enter     =  494 ! alt-enter
  integer, parameter :: curses_alt_esc       =  495 ! alt-escape
  integer, parameter :: curses_alt_bquote    =  496 ! alt-back quote
  integer, parameter :: curses_alt_lbracket  =  497 ! alt-left bracket
  integer, parameter :: curses_alt_rbracket  =  498 ! alt-right bracket
  integer, parameter :: curses_alt_semicolon =  499 ! alt-semi-colon
  integer, parameter :: curses_alt_fquote    =  500 ! alt-forward quote
  integer, parameter :: curses_alt_comma     =  501 ! alt-comma
  integer, parameter :: curses_alt_stop      =  502 ! alt-stop
  integer, parameter :: curses_alt_fslash    =  503 ! alt-forward slash
  integer, parameter :: curses_alt_bksp      =  504 ! alt-backspace
  integer, parameter :: curses_ctl_bksp      =  505 ! ctl-backspace
  integer, parameter :: curses_pad0          =  506 ! keypad 0

  integer, parameter :: curses_ctl_pad0      =  507 ! ctl-keypad 0
  integer, parameter :: curses_ctl_pad1      =  508
  integer, parameter :: curses_ctl_pad2      =  509
  integer, parameter :: curses_ctl_pad3      =  510
  integer, parameter :: curses_ctl_pad4      =  511
  integer, parameter :: curses_ctl_pad5      =  512
  integer, parameter :: curses_ctl_pad6      =  513
  integer, parameter :: curses_ctl_pad7      =  514
  integer, parameter :: curses_ctl_pad8      =  515
  integer, parameter :: curses_ctl_pad9      =  516

  integer, parameter :: curses_alt_pad0      =  517 ! alt-keypad 0
  integer, parameter :: curses_alt_pad1      =  518
  integer, parameter :: curses_alt_pad2      =  519
  integer, parameter :: curses_alt_pad3      =  520
  integer, parameter :: curses_alt_pad4      =  521
  integer, parameter :: curses_alt_pad5      =  522
  integer, parameter :: curses_alt_pad6      =  523
  integer, parameter :: curses_alt_pad7      =  524
  integer, parameter :: curses_alt_pad8      =  525
  integer, parameter :: curses_alt_pad9      =  526

  integer, parameter :: curses_ctl_del       =  527 ! clt-delete
  integer, parameter :: curses_alt_bslash    =  528 ! alt-back slash
  integer, parameter :: curses_ctl_enter     =  529 ! ctl-enter

  integer, parameter :: curses_shf_padenter  =  530 ! shift-enter on keypad
  integer, parameter :: curses_shf_padslash  =  531 ! shift-slash on keypad
  integer, parameter :: curses_shf_padstar   =  532 ! shift-star on keypad
  integer, parameter :: curses_shf_padplus   =  533 ! shift-plus on keypad
  integer, parameter :: curses_shf_padminus  =  534 ! shift-minus on keypad
  integer, parameter :: curses_shf_up        =  535 ! shift-up on keypad
  integer, parameter :: curses_shf_down      =  536 ! shift-down on keypad
  integer, parameter :: curses_shf_ic        =  537 ! shift-insert on keypad
  integer, parameter :: curses_shf_dc        =  538 ! shift-delete on keypad

  integer, parameter :: curses_key_mouse     =  539 ! "mouse" key
  integer, parameter :: curses_key_shift_l   =  540 ! Left-shift
  integer, parameter :: curses_key_shift_r   =  541 ! Right-shift
  integer, parameter :: curses_key_control_l =  542 ! Left-control
  integer, parameter :: curses_key_control_r =  543 ! Right-control
  integer, parameter :: curses_key_alt_l     =  544 ! Left-alt
  integer, parameter :: curses_key_alt_r     =  545 ! Right-alt
  integer, parameter :: curses_key_resize    =  546 ! Window resize
  integer, parameter :: curses_key_sup       =  547 ! Shifted up arrow
  integer, parameter :: curses_key_sdown     =  548 ! Shifted down arrow
  integer, parameter :: curses_key_min       =  257 ! Minimum curses key value
  integer, parameter :: curses_key_max       =  548 ! Maximum curses key value
!
!
! PDCurses Text Attributes
! ========================
!
! Originally, PDCurses used a short (16 bits) for its chtype. To include
! color, a number of things had to be sacrificed from the strict Unix and
! System V support. The main problem was fitting all character attributes
! and color into an unsigned char (all 8 bits!).
!
! Today, PDCurses by default uses a long (32 bits) for its chtype, as in
! System V. The short chtype is still available, by undefining CHTYPE_LONG
! and rebuilding the library.
!
! The following is the structure of a win->_attrs chtype:
!
! short form:
!
! -------------------------------------------------
! |15|14|13|12|11|10| 9| 8| 7| 6| 5| 4| 3| 2| 1| 0|
! -------------------------------------------------
!   color number |  attrs |   character eg 'a'
!
! The available non-color attributes are bold, reverse and blink. Others
! have no effect. The high order char is an index into an array of
! physical colors (defined in color.c) -- 32 foreground/background color
! pairs (5 bits) plus 3 bits for other attributes.
!
! long form:
!
! ----------------------------------------------------------------------------
! |31|30|29|28|27|26|25|24|23|22|21|20|19|18|17|16|15|14|13|12|..| 3| 2| 1| 0|
! ----------------------------------------------------------------------------
!       color number      |     modifiers         |      character eg 'a'
!
! The available non-color attributes are bold, underline, invisible,
! right-line, left-line, protect, reverse and blink. 256 color pairs (8
! bits), 8 bits for other attributes, and 16 bits for character data.
!
  integer(C_LONG), parameter :: curses_a_normal     = 0

  integer(C_LONG), parameter :: curses_a_altcharset = 65536
  integer(C_LONG), parameter :: curses_a_rightline  = 131072
  integer(C_LONG), parameter :: curses_a_leftline   = 262144
  integer(C_LONG), parameter :: curses_a_invis      = 524288
  integer(C_LONG), parameter :: curses_a_underline  = 1048576
  integer(C_LONG), parameter :: curses_a_reverse    = 2097152
  integer(C_LONG), parameter :: curses_a_blink      = 4194304
  integer(C_LONG), parameter :: curses_a_bold       = 8388608

  integer(C_LONG), parameter :: curses_a_attributes = -65536
  integer(C_LONG), parameter :: curses_a_chartext   = 65535
  integer(C_LONG), parameter :: curses_a_color      = -16777216

  integer(C_LONG), parameter :: curses_a_italic     = 524288
  integer(C_LONG), parameter :: curses_a_protect    = 1441792

  integer, parameter :: curses_pdc_attr_shift  = 19
  integer, parameter :: curses_pdc_color_shift = 24
!
!  VT100-compatible symbols -- box chars
!
  integer(C_LONG), parameter :: acs_ulcorner = 65644
  integer(C_LONG), parameter :: acs_llcorner = 65645
  integer(C_LONG), parameter :: acs_urcorner = 65643
  integer(C_LONG), parameter :: acs_lrcorner = 65642
  integer(C_LONG), parameter :: acs_rtee     = 65653
  integer(C_LONG), parameter :: acs_ltee     = 65652
  integer(C_LONG), parameter :: acs_btee     = 65654
  integer(C_LONG), parameter :: acs_ttee     = 65655
  integer(C_LONG), parameter :: acs_hline    = 65649
  integer(C_LONG), parameter :: acs_vline    = 65656
  integer(C_LONG), parameter :: acs_plus     = 65646
!
! VT100-compatible symbols -- other
!
  integer(C_LONG), parameter :: acs_s1       = 65647
  integer(C_LONG), parameter :: acs_s9       = 65651
  integer(C_LONG), parameter :: acs_diamond  = 65632
  integer(C_LONG), parameter :: acs_ckboard  = 65633
  integer(C_LONG), parameter :: acs_degree   = 65638
  integer(C_LONG), parameter :: acs_plminus  = 65639
  integer(C_LONG), parameter :: acs_bullet   = 65662
!
! Teletype 5410v1 symbols -- these are defined in SysV curses, but
! are not well-supported by most terminals. Stick to VT100 characters
! for optimum portability.
!
  integer(C_LONG), parameter :: acs_larrow   = 65580
  integer(C_LONG), parameter :: acs_rarrow   = 65579
  integer(C_LONG), parameter :: acs_darrow   = 65582
  integer(C_LONG), parameter :: acs_uarrow   = 65581
  integer(C_LONG), parameter :: acs_board    = 65640
  integer(C_LONG), parameter :: acs_lantern  = 65641
  integer(C_LONG), parameter :: acs_block    = 65584
!----------------------------------------------------------
!
! PDCurses External Variables
!
  integer (C_INT), bind(C, name='LINES') :: LINES
  integer (C_INT), bind(C, name='COLS') :: COLS
! windows
  type    (C_PTR), bind(C, name='stdscr') :: stdscr
  type    (C_PTR), bind(C, name='curscr') :: curscr
! screen
  type    (C_PTR), bind(C, name='sp') :: sp
  integer (C_INT), bind(C, name='COLORS') :: COLORS
  integer (C_INT), bind(C, name='COLOR_PAIRS') :: COLOR_PAIRS
  integer (C_INT), bind(C, name='TABSIZE') :: TABSIZE
! integer(C_LONG) :: acs_map
!
!----------------------------------------------------------------------
!
!   PDCurses Function Declarations
!
! Standard
!
! addch, addchnstr, addchstr, addnstr, addstr,
! attroff, attron, attrset, attr_get, attr_off,
! attr_on, attr_set,
! baudrate, beep, bkgd, bkgdset, border, box,
! can_change_color, cbreak, chgat, clearok, clear, clrtobot,
! clrtoeol, color_content, color_set, copywin, curs_set,
! def_prog_mode, def_shell_mode, delay_output,
! delch, deleteln, delscreen, delwin,
! derwin, doupdate, dupwin,
! echochar, echo, endwin, erasechar, erase,
! filter, flash, flushinp,
! getbkgd, getnstr, getstr, getwin,
! halfdelay, has_colors, has_ic, has_il, hline,
! idcok, idlok, immedok, inchnstr, inchstr, inch,
! init_color, init_pair, initscr, innstr,
! insch, insdelln, insertln, insnstr, insstr, instr,
! intrflush, isendwin, is_linetouched, is_wintouched,
! keyname, keypad, killchar, leaveok, longname,
! meta, move, mvaddch, mvaddchnstr, mvaddchstr,
! mvaddnstr, mvaddstr, mvchgat, mvcur, mvdelch,
! mvderwin, mvgetch, mvgetnstr, mvgetstr, mvhline, mvinch,
! mvinchnstr, mvinchstr, mvinnstr, mvinsch, mvinsnstr, mvinsstr,
! mvinstr, mvprintw, mvscanw, mvvline, mvwaddchnstr,
! mvwaddchstr, mvwaddch, mvwaddnstr, mvwaddstr,
! mvwchgat, mvwdelch, mvwgetch, mvwgetnstr, mvwgetstr,
! mvwhline, mvwinchnstr, mvwinchstr, mvwinch, mvwinnstr,
! mvwinsch, mvwinsnstr, mvwinsstr, mvwinstr, mvwin,
! mvwprintw, mvwscanw, mvwvline,
! napms, newpad, newterm, newwin, nl,
! nocbreak, nodelay, noecho, nonl, noqiflush, noraw, notimeout,
! overlay, overwrite,
! pair_content, pechochar, pnoutrefresh, prefresh,
! printw, putwin,
! qiflush, raw, redrawwin, refresh,
! reset_prog_mode, reset_shell_mode, resetty, ripoffline,
! savetty, scanw, scr_dump, scr_init, scr_restore, scr_set,
! scrl, scroll, scrollok, set_term, setscrreg,
! slk_attroff, slk_attr_off, slk_attron, slk_attr_on, slk_attrset,
! slk_attr_set, slk_clear, slk_color, slk_init, slk_label,
! slk_noutrefresh, slk_refresh, slk_restore, slk_set, slk_touch,
! standend, standout, start_color,
! subpad, subwin, syncok,
! termattrs, term_attrs, termname, timeout,
! touchline, touchwin, typeahead,
! untouchwin, use_env,
! vidattr, vid_attr, vidputs, vid_puts,
! vline, vw_printw, vwprintw, vw_scanw, vwscanw,
! waddchnstr, waddchstr, waddch, waddnstr, waddstr,
! wattroff, wattron, wattrset, wattr_get, wattr_off, wattr_on,
! wattr_set, wbkgdset, wbkgd, wborder, wchgat, wclear,
! wclrtobot, wclrtoeol, wcolor_set, wcursyncup, wdelch,
! wdeleteln, wechochar, werase, wgetch, wgetnstr,
! wgetstr, whline, winchnstr, winchstr, winch,
! winnstr, winsch, winsdelln, winsertln, winsnstr, winsstr,
! winstr, wmove, wnoutrefresh, wprintw, wredrawln, wrefresh,
! wscanw, wscrl, wsetscrreg, wstandend, wstandout, wsyncdown,
! wsyncup, wtimeout, wtouchln, wvline
!
  interface
!         int addch(const chtype ch);
    function addch(ch) result (ires) bind(C, name='addch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: ch
    end function addch
!         int waddch(WINDOW *win, const chtype ch);
    function waddch(win,ch) result (ires) bind(C, name='waddch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value, intent(in) :: ch
    end function waddch
!         int mvaddch(int y, int x, const chtype ch);
    function mvaddch(y,x,ch) result (ires) bind(C, name='mvaddch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value, intent(in) :: ch
    end function mvaddch
!         int mvwaddch(WINDOW *win, int y, int x, const chtype ch);
    function mvwaddch(win,y,x,ch) result (ires) bind(C, name='mvwaddch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value, intent(in) :: ch
    end function mvwaddch
!         int echochar(const chtype ch);
    function echochar(ch) result (ires) bind(C, name='echochar')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: ch
    end function echochar
!         int wechochar(WINDOW *win, const chtype ch);
    function wechochar(win,ch) result (ires) bind(C, name='wechochar')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value, intent(in) :: ch
    end function wechochar

!         int addrawch(chtype ch);
    function addrawch(ch) result (ires) bind(C, name='addrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function addrawch
!         int waddrawch(WINDOW *win, chtype ch);
    function waddrawch(win,ch) result (ires) bind(C, name='waddrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function waddrawch
!         int mvaddrawch(int y, int x, chtype ch);
    function mvaddrawch(y,x,ch) result (ires) bind(C, name='mvaddrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvaddrawch
!         int mvwaddrawch(WINDOW *win, int y, int x, chtype ch);
    function mvwaddrawch(win,y,x,ch) result (ires) bind(C, name='mvwaddrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvwaddrawch

!         int addchstr(const chtype *ch);
    function addchstr(ch) result (ires) bind(C, name='addchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), intent(in), value :: ch
    end function addchstr
!         int addchnstr(const chtype *ch, int n);
    function addchnstr(ch,n) result (ires) bind(C, name='addchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function addchnstr
!         int waddchstr(WINDOW *win, const chtype *ch);
    function waddchstr(win,ch) result (ires) bind(C, name='waddchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), intent(in), value :: ch
    end function waddchstr
!         int waddchnstr(WINDOW *win, const chtype *ch, int n);
    function waddchnstr(win,ch,n) result (ires) bind(C, name='waddchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function waddchnstr
!         int mvaddchstr(int y, int x, const chtype *ch);
    function mvaddchstr(y,x,ch) result (ires) bind(C, name='mvaddchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
    end function mvaddchstr
!         int mvaddchnstr(int y, int x, const chtype *ch, int n);
    function mvaddchnstr(y,x,ch,n) result (ires) bind(C, name='mvaddchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function mvaddchnstr
!         int mvwaddchstr(WINDOW *, int y, int x, const chtype *ch);
    function mvwaddchstr(win,y,x,ch) result (ires) bind(C, name='mvwaddchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
    end function mvwaddchstr
!         int mvwaddchnstr(WINDOW *, int y, int x, const chtype *ch, int n);
    function mvwaddchnstr(win,y,x,ch,n) result (ires) bind(C, name='mvwaddchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), intent(in), value :: ch
      integer (C_INT), value :: n
    end function mvwaddchnstr

!         int addstr(const char *str);
    function addstr(str) result (ires) bind(C, name='addstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
    end function addstr
!         int addnstr(const char *str, int n);
    function addnstr(str,n) result (ires) bind(C, name='addnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function addnstr
!         int waddstr(WINDOW *win, const char *str);
    function waddstr(win,str) result (ires) bind(C, name='waddstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
    end function waddstr
!         int waddnstr(WINDOW *win, const char *str, int n);
    function waddnstr(win,str,n) result (ires) bind(C, name='waddnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function waddnstr
!         int mvaddstr(int y, int x, const char *str);
    function mvaddstr(y,x,str) result (ires) bind(C, name='mvaddstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvaddstr
!         int mvaddnstr(int y, int x, const char *str, int n);
    function mvaddnstr(y,x,str,n) result (ires) bind(C, name='mvaddnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvaddnstr
!         int mvwaddstr(WINDOW *win, int y, int x, const char *str);
    function mvwaddstr(win,y,x,str) result (ires) bind(C, name='mvwaddstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvwaddstr
!         int mvwaddnstr(WINDOW *win, int y, int x, const char *str, int n);
    function mvwaddnstr(win,y,x,str,n) result (ires) bind(C, name='mvwaddnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvwaddnstr

!         int attroff(chtype attrs);
    function attroff(attrs) result (ires) bind(C, name='attroff')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
    end function attroff
!         int wattroff(WINDOW *win, chtype attrs);
    function wattroff(win,attrs) result (ires) bind(C, name='wattroff')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
    end function wattroff
!         int attron(chtype attrs);
    function attron(attrs) result (ires) bind(C, name='attron')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
    end function attron
!         int wattron(WINDOW *win, chtype attrs);
    function wattron(win,attrs) result (ires) bind(C, name='wattron')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
    end function wattron
!         int attrset(chtype attrs);
    function attrset(attrs) result (ires) bind(C, name='attrset')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
    end function attrset
!         int wattrset(WINDOW *win, chtype attrs);
    function wattrset(win,attrs) result (ires) bind(C, name='wattrset')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
    end function wattrset
!         int standend(void);
    function standend() result (ires) bind(C, name='standend')
      use iso_c_binding
      integer (C_INT) :: ires
    end function standend
!         int wstandend(WINDOW *win);
    function wstandend(win) result (ires) bind(C, name='wstandend')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wstandend
!         int standout(void);
    function standout() result (ires) bind(C, name='standout')
      use iso_c_binding
      integer (C_INT) :: ires
    end function standout
!         int wstandout(WINDOW *win);
    function wstandout(win) result (ires) bind(C, name='wstandout')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wstandout

!         int color_set(short color_pair, void *opts);
    function color_set(color_pair,opts) result (ires) bind(C, name='color_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function color_set
!         int wcolor_set(WINDOW *win, short color_pair, void *opts);
    function wcolor_set(win,color_pair,opts) result (ires) bind(C, name='wcolor_set')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function wcolor_set

!         int attr_get(attr_t *attrs, short *color_pair, void *opts);
    function attr_get(attrs,color_pair,opts) result (ires) bind(C, name='attr_get')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: attrs
      type (C_PTR), value :: color_pair
      type (C_PTR), value :: opts
    end function attr_get
!         int attr_off(attr_t attrs, void *opts);
    function attr_off(attrs,opts) result (ires) bind(C, name='attr_off')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function attr_off
!         int attr_on(attr_t attrs, void *opts);
    function attr_on(attrs,opts) result (ires) bind(C, name='attr_on')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function attr_on
!         int attr_set(attr_t attrs, short color_pair, void *opts);
    function attr_set(attrs,color_pair,opts) result (ires) bind(C, name='attr_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: attrs
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function attr_set
!         int wattr_get(WINDOW *win, attr_t *attrs, short *color_pair, void *opts);
    function wattr_get(win,attrs,color_pair,opts) result (ires) bind(C, name='wattr_get')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), value :: attrs
      type (C_PTR), value :: color_pair
      type (C_PTR), value :: opts
    end function wattr_get
!         int wattr_off(WINDOW *win, attr_t attrs, void *opts);
    function wattr_off(win,attrs,opts) result (ires) bind(C, name='wattr_off')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function wattr_off
!         int wattr_on(WINDOW *win, attr_t attrs, void *opts);
    function wattr_on(win,attrs,opts) result (ires) bind(C, name='wattr_on')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
      type (C_PTR), value :: opts
    end function wattr_on
!         int wattr_set(WINDOW *win, attr_t attrs, short color_pair, void *opts);
    function wattr_set(win,attrs,color_pair,opts) result (ires) bind(C, name='wattr_set')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: attrs
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function wattr_set

!         int chgat(int n, attr_t attr, short color, const void *opts);
    function chgat(n,attr,color,opts) result (ires) bind(C, name='chgat')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function chgat
!         int mvchgat(int y, int x, int n, attr_t attr, short color, const void *opts);
    function mvchgat(y,x,n,attr,color,opts) result (ires) bind(C, name='mvchgat')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function mvchgat
!         int mvwchgat(WINDOW *win, int y, int x, int n, attr_t attr, short color, const void *opts);
    function mvwchgat(win,y,x,n,attr,color,opts) result (ires) bind(C, name='mvwchgat')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function mvwchgat
!         int wchgat(WINDOW *win, int n, attr_t attr, short color, const void *opts);
    function wchgat(win,n,attr,color,opts) result (ires) bind(C, name='wchgat')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: n
      integer (C_LONG), value :: attr
      integer (C_SHORT), value :: color
      type (C_PTR), intent(in), value :: opts
    end function wchgat

!         chtype getattrs(WINDOW *win);
    function getattrs(win) result (ires) bind(C, name='getattrs')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getattrs
!         int     getbegx(WINDOW *);
    function getbegx(win) result (ires) bind(C, name='getbegx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getbegx
!         int     getbegy(WINDOW *);
    function getbegy(win) result (ires) bind(C, name='getbegy')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getbegy
!         int     getmaxx(WINDOW *);
    function getmaxx(win) result (ires) bind(C, name='getmaxx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getmaxx
!         int     getmaxy(WINDOW *);
    function getmaxy(win) result (ires) bind(C, name='getmaxy')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getmaxy
!         int     getparx(WINDOW *);
    function getparx(win) result (ires) bind(C, name='getparx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getparx
!         int     getpary(WINDOW *);
    function getpary(win) result (ires) bind(C, name='getpary')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getpary
!         int     getcurx(WINDOW *);
    function getcurx(win) result (ires) bind(C, name='getcurx')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getcurx
!         int     getcury(WINDOW *);
    function getcury(win) result (ires) bind(C, name='getcury')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getcury
!         int beep(void);
    function beep() result (ires) bind(C, name='beep')
      use iso_c_binding
      integer (C_INT) :: ires
    end function beep
!         int flash(void);
    function flash() result (ires) bind(C, name='flash')
      use iso_c_binding
      integer (C_INT) :: ires
    end function flash
!         int bkgd(chtype ch);
    function bkgd(ch) result (ires) bind(C, name='bkgd')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function bkgd
!         void bkgdset(chtype ch);
    subroutine bkgdset(ch) bind(C, name='bkgdset')
      use iso_c_binding
      integer (C_LONG), value :: ch
    end subroutine bkgdset
!         chtype getbkgd(WINDOW *win);
    function getbkgd(win) result (ires) bind(C, name='getbkgd')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function getbkgd
!         int wbkgd(WINDOW *win, chtype ch);
    function wbkgd(win,ch) result (ires) bind(C, name='wbkgd')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function wbkgd
!         void wbkgdset(WINDOW *win, chtype ch);
    subroutine wbkgdset(win,ch) bind(C, name='wbkgdset')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end subroutine wbkgdset
!         int box(WINDOW *win, chtype verch, chtype horch);
    function box(win,verch,horch) result (ires) bind(C, name='box')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: verch
      integer (C_LONG), value :: horch
    end function box
!         int hline(chtype ch, int n);
    function hline(ch,n) result (ires) bind(C, name='hline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function hline
!         int vline(chtype ch, int n);
    function vline(ch,n) result (ires) bind(C, name='vline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function vline
!         int whline(WINDOW *win, chtype ch, int n);
    function whline(win,ch,n) result (ires) bind(C, name='whline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function whline
!         int wvline(WINDOW *win, chtype ch, int n);
    function wvline(win,ch,n) result (ires) bind(C, name='wvline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function wvline
!         int mvhline(int y, int x, chtype ch, int n);
    function mvhline(y,x,ch,n) result (ires) bind(C, name='mvhline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvhline
!         int mvvline(int y, int x, chtype ch, int n);
    function mvvline(y,x,ch,n) result (ires) bind(C, name='mvvline')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvvline
!         int mvwhline(WINDOW *win, int y, int x, chtype ch, int n);
    function mvwhline(win,y,x,ch,n) result (ires) bind(C, name='mvwhline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvwhline
!         int mvwvline(WINDOW *win, int y, int x, chtype ch, int n);
    function mvwvline(win,y,x,ch,n) result (ires) bind(C, name='mvwvline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
      integer (C_INT), value :: n
    end function mvwvline
!         int clear(void);
    function clear() result (ires) bind(C, name='clear')
      use iso_c_binding
      integer (C_INT) :: ires
    end function clear
!         int wclear(WINDOW *win);
    function wclear(win) result (ires) bind(C, name='wclear')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wclear
!         int erase(void);
    function erase() result (ires) bind(C, name='erase')
      use iso_c_binding
      integer (C_INT) :: ires
    end function erase
!         int werase(WINDOW *win);
    function werase(win) result (ires) bind(C, name='werase')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function werase
!         int clrtobot(void);
    function clrtobot() result (ires) bind(C, name='clrtobot')
      use iso_c_binding
      integer (C_INT) :: ires
    end function clrtobot
!         int wclrtobot(WINDOW *win);
    function wclrtobot(win) result (ires) bind(C, name='wclrtobot')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wclrtobot
!         int clrtoeol(void);
    function clrtoeol() result (ires) bind(C, name='clrtoeol')
      use iso_c_binding
      integer (C_INT) :: ires
    end function clrtoeol
!         int wclrtoeol(WINDOW *win);
    function wclrtoeol(win) result (ires) bind(C, name='wclrtoeol')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wclrtoeol
!         int start_color(void);
    function start_color() result (ires) bind(C, name='start_color')
      use iso_c_binding
      integer (C_INT) :: ires
    end function start_color
!         int init_pair(short pair, short fg, short bg);
    function init_pair(pair,fg,bg) result (ires) bind(C, name='init_pair')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: pair
      integer (C_SHORT), value :: fg
      integer (C_SHORT), value :: bg
    end function init_pair
!         int init_color(short color, short red, short green, short blue);
    function init_color(color,red,green,blue) result (ires) bind(C, name='init_color')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color
      integer (C_SHORT), value :: red
      integer (C_SHORT), value :: green
      integer (C_SHORT), value :: blue
    end function init_color
!         bool has_colors(void);
    function has_colors() result (ires) bind(C, name='has_colors')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function has_colors
!         bool can_change_color(void);
    function can_change_color() result (ires) bind(C, name='can_change_color')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function can_change_color
!         int color_content(short color, short *red, short *green, short *blue);
    function color_content(color,red,green,blue) result (ires) bind(C, name='color_content')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color
      type (C_PTR), value :: red
      type (C_PTR), value :: green
      type (C_PTR), value :: blue
    end function color_content
!         int pair_content(short pair, short *fg, short *bg);
    function pair_content(pair,fg,bg) result (ires) bind(C, name='pair_content')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: pair
      type (C_PTR), value :: fg
      type (C_PTR), value :: bg
    end function pair_content

!         int assume_default_colors(int f, int b);
    function assume_default_colors(f,b) result (ires) bind(C, name='assume_default_colors')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: f
      integer (C_INT), value :: b
    end function assume_default_colors
!         int use_default_colors(void);
    function use_default_colors() result (ires) bind(C, name='use_default_colors')
      use iso_c_binding
      integer (C_INT) :: ires
    end function use_default_colors

!         int PDC_set_line_color(short color);
    function PDC_set_line_color(color) result (ires) bind(C, name='PDC_set_line_color')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color
    end function PDC_set_line_color
!         int delch(void);
    function delch() result (ires) bind(C, name='delch')
      use iso_c_binding
      integer (C_INT) :: ires
    end function delch
!         int wdelch(WINDOW *win);
    function wdelch(win) result (ires) bind(C, name='wdelch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wdelch
!         int mvdelch(int y, int x);
    function mvdelch(y,x) result (ires) bind(C, name='mvdelch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvdelch
!         int mvwdelch(WINDOW *win, int y, int x);
    function mvwdelch(win,y,x) result (ires) bind(C, name='mvwdelch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwdelch
!         int deleteln(void);
    function deleteln() result (ires) bind(C, name='deleteln')
      use iso_c_binding
      integer (C_INT) :: ires
    end function deleteln
!         int wdeleteln(WINDOW *win);
    function wdeleteln(win) result (ires) bind(C, name='wdeleteln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wdeleteln
!         int insdelln(int n);
    function insdelln(n) result (ires) bind(C, name='insdelln')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: n
    end function insdelln
!         int winsdelln(WINDOW *win, int n);
    function winsdelln(win,n) result (ires) bind(C, name='winsdelln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: n
    end function winsdelln
!         int insertln(void);
    function insertln() result (ires) bind(C, name='insertln')
      use iso_c_binding
      integer (C_INT) :: ires
    end function insertln
!         int winsertln(WINDOW *win);
    function winsertln(win) result (ires) bind(C, name='winsertln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function winsertln

!         int mvdeleteln(int y, int x);
    function mvdeleteln(y,x) result (ires) bind(C, name='mvdeleteln')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvdeleteln
!         int mvwdeleteln(WINDOW *win, int y, int x);
    function mvwdeleteln(win,y,x) result (ires) bind(C, name='mvwdeleteln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwdeleteln
!         int mvinsertln(int y, int x);
    function mvinsertln(y,x) result (ires) bind(C, name='mvinsertln')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvinsertln
!         int mvwinsertln(WINDOW *win, int y, int x);
    function mvwinsertln(win,y,x) result (ires) bind(C, name='mvwinsertln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwinsertln
!         int wgetch(WINDOW *win);
    function wgetch(win) result (ires) bind(C, name='wgetch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wgetch
!         int mvgetch(int y, int x);
    function mvgetch(y,x) result (ires) bind(C, name='mvgetch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvgetch
!         int mvwgetch(WINDOW *win, int y, int x);
    function mvwgetch(win,y,x) result (ires) bind(C, name='mvwgetch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwgetch
!         int ungetch(int ch);
    function ungetch(ch) result (ires) bind(C, name='ungetch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ch
    end function ungetch
!         int flushinp(void);
    function flushinp() result (ires) bind(C, name='flushinp')
      use iso_c_binding
      integer (C_INT) :: ires
    end function flushinp

!         int getstr(char *str);
    function getstr(str) result (ires) bind(C, name='getstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
    end function getstr
!         int wgetstr(WINDOW *win, char *str);
    function wgetstr(win,str) result (ires) bind(C, name='wgetstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
    end function wgetstr
!         int mvgetstr(int y, int x, char *str);
    function mvgetstr(y,x,str) result (ires) bind(C, name='mvgetstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvgetstr
!         int mvwgetstr(WINDOW *win, int y, int x, char *str);
    function mvwgetstr(win,y,x,str) result (ires) bind(C, name='mvwgetstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvwgetstr
!         int getnstr(char *str, int n);
    function getnstr(str,n) result (ires) bind(C, name='getnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function getnstr
!         int wgetnstr(WINDOW *win, char *str, int n);
    function wgetnstr(win,str,n) result (ires) bind(C, name='wgetnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function wgetnstr
!         int mvgetnstr(int y, int x, char *str, int n);
    function mvgetnstr(y,x,str,n) result (ires) bind(C, name='mvgetnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvgetnstr
!         int mvwgetnstr(WINDOW *win, int y, int x, char *str, int n);
    function mvwgetnstr(win,y,x,str,n) result (ires) bind(C, name='mvwgetnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvwgetnstr
!
! needs to be defined (is macro)
!
!         void getyx(WINDOW *win, int y, int x);
    subroutine getyx(win,y,x) bind(C, name='getyx')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end subroutine getyx

!         chtype inch(void);
    function inch() result (ires) bind(C, name='inch')
      use iso_c_binding
      integer (C_LONG) :: ires
    end function inch
!         chtype winch(WINDOW *win);
    function winch(win) result (ires) bind(C, name='winch')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
    end function winch
!         chtype mvinch(int y, int x);
    function mvinch(y,x) result (ires) bind(C, name='mvinch')
      use iso_c_binding
      integer (C_LONG) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvinch
!         chtype mvwinch(WINDOW *win, int y, int x);
    function mvwinch(win,y,x) result (ires) bind(C, name='mvwinch')
      use iso_c_binding
      integer (C_LONG) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwinch

!         int inchstr(chtype *ch);
    function inchstr(ch) result (ires) bind(C, name='inchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: ch
    end function inchstr
!         int inchnstr(chtype *ch, int n);
    function inchnstr(ch,n) result (ires) bind(C, name='inchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function inchnstr
!         int winchstr(WINDOW *win, chtype *ch);
    function winchstr(win,ch) result (ires) bind(C, name='winchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), value :: ch
    end function winchstr
!         int winchnstr(WINDOW *win, chtype *ch, int n);
    function winchnstr(win,ch,n) result (ires) bind(C, name='winchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function winchnstr
!         int mvinchstr(int y, int x, chtype *ch);
    function mvinchstr(y,x,ch) result (ires) bind(C, name='mvinchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
    end function mvinchstr
!         int mvinchnstr(int y, int x, chtype *ch, int n);
    function mvinchnstr(y,x,ch,n) result (ires) bind(C, name='mvinchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function mvinchnstr
!         int mvwinchstr(WINDOW *, int y, int x, chtype *ch);
    function mvwinchstr(win,y,x,ch) result (ires) bind(C, name='mvwinchstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
    end function mvwinchstr
!         int mvwinchnstr(WINDOW *, int y, int x, chtype *ch, int n);
    function mvwinchnstr(win,y,x,ch,n) result (ires) bind(C, name='mvwinchnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      type (C_PTR), value :: ch
      integer (C_INT), value :: n
    end function mvwinchnstr
!         WINDOW *initscr(void);
    function initscr() result (iwin) bind(C, name='initscr')
      use iso_c_binding
      type (C_PTR) :: iwin
    end function initscr
!         WINDOW *Xinitscr(int argc, char *argv[]);
!   function Xinitscr(argc,argv) result (iwin) bind(C, name='Xinitscr')
!     use iso_c_binding
!     type (C_PTR) :: iwin
!     integer (C_INT), value :: argc
!     character (C_CHAR), dimension(*) :: argv
!   end function Xinitscr
!         int endwin(void);
    function endwin() result (ires) bind(C, name='endwin')
      use iso_c_binding
      integer (C_INT) :: ires
    end function endwin
!         bool isendwin(void);
    function isendwin() result (ires) bind(C, name='isendwin')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function isendwin
!         int resize_term(int nlines, int ncols);
    function resize_term(nlines,ncols) result (ires) bind(C, name='resize_term')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
    end function resize_term
!         bool is_termresized(void);
    function is_termresized() result (ires) bind(C, name='is_termresized')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function is_termresized
!         int cbreak(void);
    function cbreak() result (ires) bind(C, name='cbreak')
      use iso_c_binding
      integer (C_INT) :: ires
    end function cbreak
!         int nocbreak(void);
    function nocbreak() result (ires) bind(C, name='nocbreak')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nocbreak
!         int echo(void);
    function echo() result (ires) bind(C, name='echo')
      use iso_c_binding
      integer (C_INT) :: ires
    end function echo
!         int noecho(void);
    function noecho() result (ires) bind(C, name='noecho')
      use iso_c_binding
      integer (C_INT) :: ires
    end function noecho
!         int halfdelay(int tenths);
    function halfdelay(tenths) result (ires) bind(C, name='halfdelay')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: tenths
    end function halfdelay
!         int intrflush(WINDOW *win, bool bf);
    function intrflush(win, bf) result (ires) bind(C, name='intrflush')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function intrflush
!         int keypad(WINDOW *win, bool bf);
    function keypad(win, bf) result (ires) bind(C, name='keypad')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function keypad
!         int meta(WINDOW *win, bool bf);
    function meta(win, bf) result (ires) bind(C, name='meta')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function meta
!         int nl(void);
    function nl() result (ires) bind(C, name='nl')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nl
!         int nonl(void);
    function nonl() result (ires) bind(C, name='nonl')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nonl
!         int nodelay(WINDOW *win, bool bf);
    function nodelay(win, bf) result (ires) bind(C, name='nodelay')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function nodelay
!         int notimeout(WINDOW *win, bool bf);
    function notimeout(win, bf) result (ires) bind(C, name='notimeout')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function notimeout
!         int raw(void);
    function raw() result (ires) bind(C, name='raw')
      use iso_c_binding
      integer (C_INT) :: ires
    end function raw
!         int noraw(void);
    function noraw() result (ires) bind(C, name='noraw')
      use iso_c_binding
      integer (C_INT) :: ires
    end function noraw
!         void noqiflush(void);
    subroutine noqiflush() bind(C, name='noqiflush')
      use iso_c_binding
    end subroutine noqiflush
!         void qiflush(void);
    subroutine qiflush() bind(C, name='qiflush')
      use iso_c_binding
    end subroutine qiflush
!         void timeout(int delay);
    subroutine timeout(delay) bind(C, name='timeout')
      use iso_c_binding
      integer (C_INT), value :: delay
    end subroutine timeout
!         void wtimeout(WINDOW *win, int delay);
    subroutine wtimeout(win,delay) bind(C, name='wtimeout')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_INT), value :: delay
    end subroutine wtimeout
!         int typeahead(int fildes);
    function typeahead(fildes) result (ires) bind(C, name='typeahead')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: fildes
    end function typeahead

!         int crmode(void);
    function crmode() result (ires) bind(C, name='crmode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function crmode
!         int nocrmode(void);
    function nocrmode() result (ires) bind(C, name='nocrmode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function nocrmode
!         int insch(chtype ch);
    function insch(ch) result (ires) bind(C, name='insch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function insch
!         int winsch(WINDOW *win, chtype ch);
    function winsch(win,ch) result (ires) bind(C, name='winsch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function winsch
!         int mvinsch(int y, int x, chtype ch);
    function mvinsch(y,x,ch) result (ires) bind(C, name='mvinsch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvinsch
!         int mvwinsch(WINDOW *win, int y, int x, chtype ch);
    function mvwinsch(win,y,x,ch) result (ires) bind(C, name='mvwinsch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvwinsch

!         int insrawch(chtype ch);
    function insrawch(ch) result (ires) bind(C, name='insrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value :: ch
    end function insrawch
!         int winsrawch(WINDOW *win, chtype ch);
    function winsrawch(win,ch) result (ires) bind(C, name='winsrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_LONG), value :: ch
    end function winsrawch
!         int mvinsrawch(int y, int x, chtype ch);
    function mvinsrawch(y,x,ch) result (ires) bind(C, name='mvinsrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvinsrawch
!         int mvwinsrawch(WINDOW *win, int y, int x, chtype ch);
    function mvwinsrawch(win,y,x,ch) result (ires) bind(C, name='mvwinsrawch')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      integer (C_LONG), value :: ch
    end function mvwinsrawch
!         int insstr(const char *str);
    function insstr(str) result (ires) bind(C, name='insstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
    end function insstr
!         int insnstr(const char *str, int n);
    function insnstr(str,n) result (ires) bind(C, name='insnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function insnstr
!         int winsstr(WINDOW *win, const char *str);
    function winsstr(win,str) result (ires) bind(C, name='winsstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
    end function winsstr
!         int winsnstr(WINDOW *win, const char *str, int n);
    function winsnstr(win,str,n) result (ires) bind(C, name='winsnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function winsnstr
!         int mvinsstr(int y, int x, const char *str);
    function mvinsstr(y,x,str) result (ires) bind(C, name='mvinsstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvinsstr
!         int mvinsnstr(int y, int x, const char *str, int n);
    function mvinsnstr(y,x,str,n) result (ires) bind(C, name='mvinsnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvinsnstr
!         int mvwinsstr(WINDOW *win, int y, int x, const char *str);
    function mvwinsstr(win,y,x,str) result (ires) bind(C, name='mvwinsstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
    end function mvwinsstr
!         int mvwinsnstr(WINDOW *win, int y, int x, const char *str, int n);
    function mvwinsnstr(win,y,x,str,n) result (ires) bind(C, name='mvwinsnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: str
      integer (C_INT), value :: n
    end function mvwinsnstr

!         int instr(char *str);
    function instr(str) result (ires) bind(C, name='instr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
    end function instr
!         int innstr(char *str, int n);
    function innstr(str,n) result (ires) bind(C, name='innstr')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function innstr
!         int winstr(WINDOW *win, char *str);
    function winstr(win,str) result (ires) bind(C, name='winstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
    end function winstr
!         int winnstr(WINDOW *win, char *str, int n);
    function winnstr(win,str,n) result (ires) bind(C, name='winnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function winnstr
!         int mvinstr(int y, int x, char *str);
    function mvinstr(y,x,str) result (ires) bind(C, name='mvinstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvinstr
!         int mvinnstr(int y, int x, char *str, int n);
    function mvinnstr(y,x,str,n) result (ires) bind(C, name='mvinnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvinnstr
!         int mvwinstr(WINDOW *win, int y, int x, char *str);
    function mvwinstr(win,y,x,str) result (ires) bind(C, name='mvwinstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
    end function mvwinstr
!         int mvwinnstr(WINDOW *win, int y, int x, char *str, int n);
    function mvwinnstr(win,y,x,str,n) result (ires) bind(C, name='mvwinnstr')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*) :: str
      integer (C_INT), value :: n
    end function mvwinnstr
!         int def_prog_mode(void);
    function def_prog_mode() result (ires) bind(C, name='def_prog_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function def_prog_mode
!         int def_shell_mode(void);
    function def_shell_mode() result (ires) bind(C, name='def_shell_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function def_shell_mode
!         int reset_prog_mode(void);
    function reset_prog_mode() result (ires) bind(C, name='reset_prog_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function reset_prog_mode
!         int reset_shell_mode(void);
    function reset_shell_mode() result (ires) bind(C, name='reset_shell_mode')
      use iso_c_binding
      integer (C_INT) :: ires
    end function reset_shell_mode
!         int resetty(void);
    function resetty() result (ires) bind(C, name='resetty')
      use iso_c_binding
      integer (C_INT) :: ires
    end function resetty
!         int savetty(void);
    function savetty() result (ires) bind(C, name='savetty')
      use iso_c_binding
      integer (C_INT) :: ires
    end function savetty
!         int ripoffline(int line, int (*init)(WINDOW *, int));
!         int copywin(const WINDOW *, WINDOW *, int, int, int, int, int, int, int);
    function copywin(srcwin, dstwin, sminrow, smincol, dminrow, dmincol,  &
                     dmaxrow, dmaxcol, overlay) result (ires) bind(C, name='copywin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value, intent(in) :: srcwin
      type (C_PTR), value :: dstwin
      integer (C_INT), value :: sminrow, smincol
      integer (C_INT), value :: dminrow, dmincol
      integer (C_INT), value :: dmaxrow, dmaxcol
      integer (C_INT), value :: overlay
    end function copywin
!         int overlay(const WINDOW *srcwin, WINDOW *dstwin);
    function overlay(srcwin, dstwin) result (ires) bind(C, name='overlay')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value, intent(in) :: srcwin
      type (C_PTR), value :: dstwin
    end function overlay
!         int overwrite(const WINDOW *srcwin, WINDOW *dstwin);
    function overwrite(srcwin, dstwin) result (ires) bind(C, name='overwrite')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value, intent(in) :: srcwin
      type (C_PTR), value :: dstwin
    end function overwrite
!         int curs_set(int visibility);
    function curs_set(visibility) result (ires) bind(C, name='curs_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: visibility
    end function curs_set
!         int napms(int ms);
    function napms(ms) result (ires) bind(C, name='napms')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ms
    end function napms

!         int draino(int ms);
    function draino(ms) result (ires) bind(C, name='draino')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ms
    end function draino
!         int resetterm(void);
    function resetterm() result (ires) bind(C, name='resetterm')
      use iso_c_binding
      integer (C_INT) :: ires
    end function resetterm
!         int fixterm(void);
    function fixterm() result (ires) bind(C, name='fixterm')
      use iso_c_binding
      integer (C_INT) :: ires
    end function fixterm
!         int saveterm(void);
    function saveterm() result (ires) bind(C, name='saveterm')
      use iso_c_binding
      integer (C_INT) :: ires
    end function saveterm
!         char *keyname(int key);
    function keyname(key) result (ich) bind(C, name='keyname')
      use iso_c_binding
      character (C_CHAR) :: ich
      integer (C_INT), value :: key
    end function keyname
!         bool has_key(int key);
    function has_key(key) result (ires) bind(C, name='has_key')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      integer (C_INT), value :: key
    end function has_key
!         int mouse_set(unsigned long mbe);
    function mouse_set(mbe) result (ires) bind(C, name='mouse_set')
      use iso_c_binding
      integer (C_LONG) :: mbe
      integer (C_INT) :: ires
    end function mouse_set
!         int mouse_on(unsigned long mbe);
    function mouse_on(mbe) result (ires) bind(C, name='mouse_on')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG) :: mbe
    end function mouse_on
!         int mouse_off(unsigned long mbe);
    function mouse_off(mbe) result (ires) bind(C, name='mouse_off')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG) :: mbe
    end function mouse_off
!         int request_mouse_pos(void);
    function request_mouse_pos() result (ires) bind(C, name='request_mouse_pos')
      use iso_c_binding
      integer (C_INT) :: ires
    end function request_mouse_pos
!         int map_button(unsigned long button);
    function map_button(button) result (ires) bind(C, name='map_button')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG) :: button
    end function map_button
!         void wmouse_position(WINDOW *win, int *y, int *x);
    subroutine wmouse_position(win,y,x) bind(C, name='wmouse_position')
      use iso_c_binding
      type (C_PTR), value :: win
      type (C_PTR), value :: y
      type (C_PTR), value :: x
    end subroutine wmouse_position
!         unsigned long getmouse(void);
!         bool wenclose(const WINDOW *win, int y, int x);
    function wenclose(win,y,x) result (ires) bind(C, name='wenclose')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      type (C_PTR), intent(in) :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function wenclose
!         int move(int y, int x);
    function move(y,x) result (ires) bind(C, name='move')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function move
!         int wmove(WINDOW *win, int y, int x);
    function wmove(win,y,x) result (ires) bind(C, name='wmove')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function wmove

!         int clearok(WINDOW *win, bool bf);
    function clearok(win, bf) result (ires) bind(C, name='clearok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function clearok
!         int idlok(WINDOW *win, bool bf);
    function idlok(win, bf) result (ires) bind(C, name='idlok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function idlok
!         void idcok(WINDOW *win, bool bf);
    subroutine idcok(win, bf) bind(C, name='idcok')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end subroutine idcok
!         void immedok(WINDOW *win, bool bf);
    subroutine immedok(win, bf) bind(C, name='immedok')
      use iso_c_binding
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end subroutine immedok
!         int leaveok(WINDOW *win, bool bf);
    function leaveok(win, bf) result (ires) bind(C, name='leaveok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function leaveok
!         int setscrreg(int top, int bot);
    function setscrreg(top,bot) result (ires) bind(C, name='setscrreg')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: top
      integer (C_INT), value :: bot
    end function setscrreg
!         int wsetscrreg(WINDOW *win, int top, int bot);
    function wsetscrreg(win,top,bot) result (ires) bind(C, name='wsetscrreg')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: top
      integer (C_INT), value :: bot
    end function wsetscrreg
!         int scrollok(WINDOW *win, bool bf);
    function scrollok(win, bf) result (ires) bind(C, name='scrollok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function scrollok

!         int raw_output(bool bf);
    function raw_output(bf) result (ires) bind(C, name='raw_output')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SIGNED_CHAR), value :: bf
    end function raw_output
!         int refresh(void);
    function refresh() result (ires) bind(C, name='refresh')
      use iso_c_binding
      integer (C_INT) :: ires
    end function refresh
!         int wrefresh(WINDOW *win);
    function wrefresh(win) result (ires) bind(C, name='wrefresh')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wrefresh
!         int wnoutrefresh(WINDOW *win);
    function wnoutrefresh(win) result (ires) bind(C, name='wnoutrefresh')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function wnoutrefresh
!         int doupdate(void);
    function doupdate() result (ires) bind(C, name='doupdate')
      use iso_c_binding
      integer (C_INT) :: ires
    end function doupdate
!         int redrawwin(WINDOW *win);
    function redrawwin(win) result (ires) bind(C, name='redrawwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function redrawwin
!         int wredrawln(WINDOW *win, int beg_line, int num_lines);
    function wredrawln(win,beg_line,num_lines) result (ires) bind(C, name='wredrawln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: beg_line
      integer (C_INT), value :: num_lines
    end function wredrawln

!         int scanw(const char *fmt, ...);
    function scanw(formt) result (ires) bind(C, name='scanw')
      use iso_c_binding
      integer (C_INT) :: ires
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function scanw
!         int wscanw(WINDOW *win, const char *fmt, ...);
    function wscanw(win,formt) result (ires) bind(C, name='wscanw')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function wscanw
!         int mvscanw(int y, int x, const char *fmt, ...);
    function mvscanw(y,x,formt) result (ires) bind(C, name='mvscanw')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function mvscanw
!         int mvwscanw(WINDOW *win, int y, int x, const char *fmt, ...);
    function mvwscanw(win,y,x,formt) result (ires) bind(C, name='mvwscanw')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
      character (C_CHAR), dimension(*), intent(in) :: formt
    end function mvwscanw
!         int vwscanw(WINDOW *win, const char *fmt, va_list varglist);
!         int vw_scanw(WINDOW *win, const char *fmt, va_list varglist);
!         int scroll(WINDOW *win);
    function scroll(win) result (ires) bind(C, name='scroll')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function scroll
!         int scrl(int n);
    function scrl(n) result (ires) bind(C, name='scrl')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: n
    end function scrl
!         int wscrl(WINDOW *win, int n);
    function wscrl(win,n) result (ires) bind(C, name='wscrl')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: n
    end function wscrl

!         void filter(void);
    subroutine filter() bind(C, name='filter')
      use iso_c_binding
    end subroutine filter
!         void use_env(bool x);
    subroutine use_env(x) bind(C, name='use_env')
      use iso_c_binding
      integer (C_SIGNED_CHAR), value :: x
    end subroutine use_env
!         int delay_output(int ms);
    function delay_output(ms) result (ires) bind(C, name='delay_output')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: ms
    end function delay_output
!         int slk_init(int fmt);
    function slk_init(fmt) result (ires) bind(C, name='slk_init')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: fmt
    end function slk_init
!         int slk_set(int labnum, const char *label, int justify);
    function slk_set(labnum,label,justify) result (ires) bind(C, name='slk_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: labnum
      character (C_CHAR), dimension(*), intent(in) :: label
      integer (C_INT), value :: justify
    end function slk_set
!         int slk_refresh(void);
    function slk_refresh() result (ires) bind(C, name='slk_refresh')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_refresh
!         int slk_noutrefresh(void);
    function slk_noutrefresh() result (ires) bind(C, name='slk_noutrefresh')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_noutrefresh
!         char *slk_label(int labnum);
!         int slk_clear(void);
    function slk_clear() result (ires) bind(C, name='slk_clear')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_clear
!         int slk_restore(void);
    function slk_restore() result (ires) bind(C, name='slk_restore')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_restore
!         int slk_touch(void);
    function slk_touch() result (ires) bind(C, name='slk_touch')
      use iso_c_binding
      integer (C_INT) :: ires
    end function slk_touch
!         int slk_attron(const chtype attrs);
    function slk_attron(attrs) result (ires) bind(C, name='slk_attron')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
    end function slk_attron
!         int slk_attr_on(const attr_t attrs, void *opts);
    function slk_attr_on(attrs,opts) result (ires) bind(C, name='slk_attr_on')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
      type (C_PTR), value :: opts
    end function slk_attr_on
!         int slk_attrset(const chtype attrs);
    function slk_attrset(attrs) result (ires) bind(C, name='slk_attrset')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
    end function slk_attrset
!         int slk_attr_set(const attr_t attrs, short color_pair, void *opts);
    function slk_attr_set(attrs,color_pair,opts) result (ires) bind(C, name='slk_attr_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function slk_attr_set
!         int slk_attroff(const chtype attrs);
    function slk_attroff(attrs) result (ires) bind(C, name='slk_attroff')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
    end function slk_attroff
!         int slk_attr_off(const attr_t attrs, void *opts);
    function slk_attr_off(attrs,opts) result (ires) bind(C, name='slk_attr_off')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
      type (C_PTR), value :: opts
    end function slk_attr_off
!         int slk_color(short color_pair);
    function slk_color(color_pair) result (ires) bind(C, name='slk_color')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color_pair
    end function slk_color
!         int PDC_mouse_in_slk(int y, int x);
    function PDC_mouse_in_slk(y,x) result (ires) bind(C, name='PDC_mouse_in_slk')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function PDC_mouse_in_slk
!         void PDC_slk_free(void);
    subroutine PDC_slk_free() bind(C, name='PDC_slk_free')
      use iso_c_binding
    end subroutine PDC_slk_free
!         void PDC_slk_initialize(void);
    subroutine PDC_slk_initialize() bind(C, name='PDC_slk_initialize')
      use iso_c_binding
    end subroutine PDC_slk_initialize
!         int baudrate(void);
    function baudrate() result (ires) bind(C, name='baudrate')
      use iso_c_binding
      integer (C_INT) :: ires
    end function baudrate
!         char erasechar(void);
    function erasechar() result (ich) bind(C, name='erasechar')
      use iso_c_binding
      character (C_CHAR) :: ich
    end function erasechar
!         bool has_ic(void);
    function has_ic() result (ires) bind(C, name='has_ic')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function has_ic
!         bool has_il(void);
    function has_il() result (ires) bind(C, name='has_il')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function has_il
!         char killchar(void);
    function killchar() result (ich) bind(C, name='killchar')
      use iso_c_binding
      character (C_CHAR) :: ich
    end function killchar
!         char *longname(void);
!         chtype termattrs(void);
    function termattrs() result (ires) bind(C, name='termattrs')
      use iso_c_binding
      integer (C_LONG) :: ires
    end function termattrs
!         attr_t term_attrs(void);
    function term_attrs() result (ires) bind(C, name='term_attrs')
      use iso_c_binding
      integer (C_LONG) :: ires
    end function term_attrs
!         char *termname(void);
!         char wordchar(void);
    function wordchar() result (ich) bind(C, name='wordchar')
      use iso_c_binding
      character (C_CHAR) :: ich
    end function wordchar
!         int mvcur(int oldrow, int oldcol, int newrow, int newcol);
    function mvcur(oldrow,oldcol,newrow,newcol) result (ires) bind(C, name='mvcur')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: oldrow
      integer (C_INT), value :: oldcol
      integer (C_INT), value :: newrow
      integer (C_INT), value :: newcol
    end function mvcur

!         int touchwin(WINDOW *win);
    function touchwin(win) result (ires) bind(C, name='touchwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function touchwin
!         int touchline(WINDOW *win, int start, int count);
    function touchline(win,start,cnt) result (ires) bind(C, name='touchline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: start
      integer (C_INT), value :: cnt
    end function touchline
!         int untouchwin(WINDOW *win);
    function untouchwin(win) result (ires) bind(C, name='untouchwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function untouchwin
!         int wtouchln(WINDOW *win, int y, int n, int changed);
    function wtouchln(win,y,n,changed) result (ires) bind(C, name='wtouchln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: n
      integer (C_INT), value :: changed
    end function wtouchln
!         bool is_linetouched(WINDOW *win, int line);
    function is_linetouched(win,line) result (ires) bind(C, name='is_linetouched')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: line
    end function is_linetouched
!         bool is_wintouched(WINDOW *win);
    function is_wintouched(win) result (ires) bind(C, name='is_wintouched')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      type (C_PTR), value :: win
    end function is_wintouched
!         WINDOW *newwin(int nlines, int ncols, int begy, int begx);
    function newwin(nlines,ncols,begy,begx) result (iwin) bind(C, name='newwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function newwin
!       WINDOW *derwin(WINDOW* orig, int nlines, int ncols, int begy, int begx);
    function derwin(orig,nlines,ncols,begy,begx) result (iwin) bind(C, name='derwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: orig
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function derwin
!       WINDOW *subwin(WINDOW* orig, int nlines, int ncols, int begy, int begx);
    function subwin(orig,nlines,ncols,begy,begx) result (iwin) bind(C, name='subwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: orig
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function subwin
!         WINDOW *dupwin(WINDOW *win);
    function dupwin(win) result (iwin) bind(C, name='dupwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: win
    end function dupwin
!         int delwin(WINDOW *win);
    function delwin(win) result (ires) bind(C, name='delwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function delwin
!         int mvwin(WINDOW *win, int y, int x);
    function mvwin(win,y,x) result (ires) bind(C, name='mvwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwin
!         int mvderwin(WINDOW *win, int pary, int parx);
    function mvderwin(win,pary,parx) result (ires) bind(C, name='mvderwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: pary
      integer (C_INT), value :: parx
    end function mvderwin
!         int syncok(WINDOW *win, bool bf);
    function syncok(win, bf) result (ires) bind(C, name='syncok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function syncok
!         void wsyncup(WINDOW *win);
    subroutine wsyncup(win) bind(C, name='wsyncup')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine wsyncup
!         void wcursyncup(WINDOW *win);
    subroutine wcursyncup(win) bind(C, name='wcursyncup')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine wcursyncup
!         void wsyncdown(WINDOW *win);
    subroutine wsyncdown(win) bind(C, name='wsyncdown')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine wsyncdown

!         WINDOW *resize_window(WINDOW *win, int nlines, int ncols);
    function resize_window(win,nlines,ncols) result (iwin) bind(C, name='resize_window')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: win
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
    end function resize_window
!         int wresize(WINDOW *win, int nlines, int ncols);
    function wresize(win,nlines,ncols) result (ires) bind(C, name='wresize')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
    end function wresize
!         WINDOW *PDC_makelines(WINDOW *win);
    function PDC_makelines(win) result (iwin) bind(C, name='PDC_makelines')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: win
    end function PDC_makelines
!         WINDOW *PDC_makenew(int nlines, int ncols, int begy, int begx);
    function PDC_makenew(nlines,ncols,begy,begx) result (iwin) bind(C, name='PDC_makenew')
      use iso_c_binding
      type (C_PTR) :: iwin
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function PDC_makenew
!         void PDC_sync(WINDOW *win);
    subroutine PDC_sync(win) bind(C, name='PDC_sync')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine PDC_sync
  end interface
contains
!         int getch(void);
  function getch() result(ires)
    integer (C_INT) :: ires
    ires=wgetch(stdscr)
  end function getch
end module curses
!
! Test program, after testcurses.c
!
module winsize
  use curses
  integer (C_INT) :: height, width
end module winsize

module commands
  use curses
  integer, parameter :: MAX_OPTIONS = 8
  character (len=14), dimension(MAX_OPTIONS) :: command =  &
  (/ 'Intro test    ',  &
     'Resize test   ',  &
     'Scroll test   ',  &
     'Input test    ',  &
     'Output test   ',  &
     'ACS test      ',  &
     'Colour test   ',  &
     'Clipboard test' /)
end module commands

subroutine initTest(iwin, istat)
  use curses
  use winsize
  type (C_PTR), intent(inout) :: iwin
  integer, intent(out) :: istat
  integer (C_INT) :: ierr

  istat=0
  iwin = C_NULL_PTR
  iwin = initscr()
  ierr = start_color()
  width=60
  height=13
  iwin=newwin(height, width, (LINES-height)/2, (COLS-width)/2)
  if (.not.c_associated(iwin)) then
    istat=1
  end if
end subroutine initTest

subroutine introtest(iwin, istat)
  use curses
  use winsize
  type (C_PTR), intent(inout) :: iwin
  integer, intent(out) :: istat
  integer (C_INT) :: ierr

  istat=0
  ierr=werase(iwin)
  ierr=wmove(iwin, height/2-5, width/2-10)
  ierr=wvline(iwin, acs_vline, 10)
  ierr=wmove(iwin, height/2, width/2-10)
  ierr=whline(iwin, acs_hline, 20)
  call cont(iwin)
  ierr=beep()
  ierr=werase(iwin)
  ierr=box(iwin, acs_vline, acs_hline)
  ierr=cbreak()
  ierr=mvwaddstr(iwin, 1, 1,  &
        'You should have a rectangle in the middle of the screen' // C_NULL_CHAR)
  ierr=mvwaddstr(iwin, 2, 1, 'You should have heard a beep' // C_NULL_CHAR)
  call cont(iwin)
  ierr=flash()
  ierr=mvwaddstr(iwin, 3, 1, 'You should have seen a flash' // C_NULL_CHAR)
  call cont(iwin)
end subroutine introtest

subroutine resizetest(istat)
  use curses
  integer, intent(out) :: istat

  type (C_PTR) :: win1
  integer (C_INT) :: nwidth = 135, nheight = 52
  integer (C_INT) :: ierr, owidth, oheight
  character (len=80) :: str

  istat=0
  owidth = COLS
  oheight = LINES

  ierr=savetty()
  ierr=resize_term(nheight, nwidth)
  if (ierr /= 0) then
    write(*,*) 'ERROR: failure in resize_term'
  end if
  ierr=clear()
  ierr=refresh()

  win1 = newwin(10, 50, 14, 25)
  if (.not.c_associated(win1)) then
    ierr=endwin()
    istat=-1
    return
  end if
  ierr=mvwaddstr(win1, 0, 0, 'The screen may now be resized' // C_NULL_CHAR)
  write(str, '(a,1x,i0,1x,a,1x,i0)') 'Given size:', nwidth, 'by', nheight
  ierr=mvwaddstr(win1, 1, 4, trim(str) // C_NULL_CHAR)
  write(str, '(a,1x,i0,1x,a,1x,i0)') 'Actual size:', COLS, 'by', LINES
  ierr=mvwaddstr(win1, 1, 4, trim(str) // C_NULL_CHAR)
  call cont(win1)

  ierr=wclear(win1)
  ierr=resetty()

  ierr=mvwaddstr(win1, 0, 0, 'The screen should now be reset' // C_NULL_CHAR)
  write(str, '(a,1x,i0,1x,a,1x,i0)') 'Old size:', owidth, 'by', oheight
  ierr=mvwaddstr(win1, 1, 4, trim(str) // C_NULL_CHAR)
  write(str, '(a,1x,i0,1x,a,1x,i0)') 'Size now:', COLS, 'by', LINES
  ierr=mvwaddstr(win1, 1, 4, trim(str) // C_NULL_CHAR)
  call cont(win1)

  ierr=delwin(win1)

  ierr=clear()
  ierr=refresh()
end subroutine resizetest

subroutine scrolltest(win, istat)
  use curses
  use winsize
  type (C_PTR), intent(inout) :: win
  integer, intent(out) :: istat
  integer (C_INT) :: ierr
  integer :: i, oldy

  istat=0
  ierr=werase(win)
  ierr=mvwaddstr(win, height - 2, 1,  &
                 'The window will now scroll slowly' // C_NULL_CHAR)
  ierr=box(win, acs_vline, acs_hline)
  ierr=wrefresh(win)
  ierr=scrollok(win, curses_true)
  ierr=napms(500)

  do i = 1, height
    ierr=napms(150)
    ierr=scroll(win)
    ierr=wrefresh(win)
  end do

  oldy = getmaxy(win)
  ierr=mvwaddstr(win, 6, 1, 'The top of the window will scroll' // C_NULL_CHAR)
  ierr=wmove(win, 1, 1)
  ierr=wsetscrreg(win, 0, 4)
  ierr=box(win, acs_vline, acs_hline)
  ierr=wrefresh(win)

  do i = 1, 5
    ierr=napms(500)
    ierr=scroll(win)
    ierr=wrefresh(win)
  end do

  ierr=mvwaddstr(win, 3, 1,  &
                 'The bottom of the window will scroll' // C_NULL_CHAR)
  ierr=wmove(win, 8, 1)
  oldy=old-1
  ierr=wsetscrreg(win, 5, oldy)
  ierr=box(win, acs_vline, acs_hline)
  ierr=wrefresh(win)

  do i = 5, oldy
    ierr=napms(300)
    ierr=wscrl(win, -1)
    ierr=wrefresh(win)
  end do
  ierr=wsetscrreg(win, 0, oldy)
end subroutine scrolltest

subroutine outputtest(win)
  use curses
  type (C_PTR), intent(inout) :: win
  type (C_PTR) :: win1
  character (len=80) :: buffer
  integer (C_LONG) :: ch
  integer (C_INT) :: by, bx, ierr

  ierr=nl()
  ierr=wclear(win)
  ierr=mvwaddstr(win, 1, 1,  &
                 'You should now have a screen in the upper '  //  &
                 'left corner, and this text should have wrapped' // C_NULL_CHAR)
  ierr=waddstr(win,'\nThis text should be down\n' // C_NULL_CHAR)
  ierr=waddstr(win, 'and broken into two here ^' // C_NULL_CHAR)
  call cont(win)

  ierr=wclear(win)
  ierr=wattron(win, curses_a_bold)
  ierr=mvwaddstr(win, 1, 1, 'A new window will appear with this text in it' // C_NULL_CHAR)
  ierr=mvwaddstr(win, 8, 1, 'Press any key to continue' // C_NULL_CHAR)
  ierr=wrefresh(win)
  ierr=wgetch(win)

  if (LINES < 24 .or. COLS < 75) then
    ierr=mvwaddstr(win, 5, 1, 'Some tests have been skipped as they require a' // C_NULL_CHAR)
    ierr=mvwaddstr(win, 6, 1, 'display of at least 24 LINES by 75 COLUMNS' // C_NULL_CHAR)
    call cont(win)
  else
    win1 = newwin(10, 50, 14, 25)
    if (.not.c_associated(win1)) then
      istat=1
      ierr=endwin()
      return
    end if
    ierr=wbkgd(win1, curses_a_normal)

    ierr=wclear(win1)
    ierr=mvwaddstr(win1, 5, 1,  &
           'This text should appear; using overlay option' // C_NULL_CHAR)
    ierr=copywin(win, win1, 0, 0, 0, 0, 9, 49, 1)
    ierr=box(win1, acs_vline, acs_hline)
    ierr=wmove(win1, 8, 26)
    ierr=wrefresh(win1)
    ierr=wgetch(win1)

    ierr=wclear(win1)

    ierr=wattron(win1, curses_a_blink)
    ierr=mvwaddstr(win1, 4, 1,  &
           'This blinking text should appear in only the second window' // C_NULL_CHAR)
    ierr=wattroff(win1, curses_a_blink)

    ierr=mvwin(win1, by, bx)
    ierr=overlay(win, win1)
    ierr=mvwin(win1, 14, 25)
    ierr=wmove(win1, 8, 26)
    ierr=wrefresh(win1)
    ierr=wgetch(win1)

    ierr=delwin(win1)
  end if

  ierr=clear()
  ierr=wclear(win)
  ierr=wrefresh(win)
  ierr=mvwaddstr(win, 6, 2, 'This line shouldn''t appear' // C_NULL_CHAR)
  ierr=mvwaddstr(win, 4, 2, 'Only half of the next line is visible' // C_NULL_CHAR)
  ierr=mvwaddstr(win, 5, 2, 'Only half of the next line is visible' // C_NULL_CHAR)
  ierr=wmove(win, 6, 1)
  ierr=wclrtobot(win)
  ierr=wmove(win, 5, 20)
  ierr=wclrtoeol(win)
  ierr=mvwaddstr(win, 8, 2, 'This line also shouldn''t appear' // C_NULL_CHAR)
  ierr=wmove(win, 8, 1)
  ierr=winsdelln(win, -1)
  call cont(win)

  ierr=wmove(win, 5, 9)
  ch = winch(win)

  ierr=wclear(win)
  ierr=wmove(win, 6, 2)
  ierr=waddstr(win, 'The next char should be l:  ' // C_NULL_CHAR)
  ierr=winsch(win, ch)
  call cont(win)

  ierr=mvwinsstr(win, 6, 2, 'A1B2C3D4E5' // C_NULL_CHAR)
  call cont(win)

  ierr=wmove(win, 5, 1)
  ierr=winsdelln(win, 1)
  ierr=mvwaddstr(win, 5, 2, 'The lines below should have moved down' // C_NULL_CHAR)
  call cont(win)

  ierr=wclear(win)
  ierr=wmove(win, 2, 2)
  write(buffer, '(a,1x,i0,1x,a)')  &
    'This is a formatted string in a window:', 42, 'is it'
  ierr=mvwaddstr(win, 1, 1, trim(buffer) // C_NULL_CHAR)
  ierr=mvwaddstr(win, 10, 1, 'Enter a string: ' // C_NULL_CHAR)
  ierr=wrefresh(win)
  ierr=echo()
! ierr=wscanw(win, "%s", buffer)
! write(*,*) 'Entered: ', trim(buffer)
  ierr=mvaddstr(10, 1, 'Enter a string: ' // C_NULL_CHAR)
! ierr=scanw("%s", buffer)
! write(*,*) 'Entered: ', trim(buffer)

  ierr=wclear(win)
  ierr=curs_set(2)
  ierr=mvwaddstr(win, 1, 1, 'The cursor should be in high-visibility mode' // C_NULL_CHAR)
  call cont(win)

  ierr=wclear(win)
  ierr=curs_set(0)
  ierr=mvwaddstr(win, 1, 1, 'The cursor should have disappeared' // C_NULL_CHAR)
  call cont(win)

  ierr=wclear(win)
  ierr=curs_set(1);
  ierr=mvwaddstr(win, 1, 1, 'The cursor should be normal' // C_NULL_CHAR)
  call cont(win)
end subroutine outputtest

subroutine acstest(win)
  use curses
  type (C_PTR), intent(inout) :: win
  integer, parameter :: ACSNUM = 25
  character (len=13), dimension(ACSNUM) :: acs_names = (/  &
    'ACS_ULCORNER ', 'ACS_URCORNER ', 'ACS_LLCORNER ', 'ACS_LRCORNER ',  &
    'ACS_LTEE     ', 'ACS_RTEE     ', 'ACS_TTEE     ', 'ACS_BTEE     ',  &
    'ACS_HLINE    ', 'ACS_VLINE    ', 'ACS_PLUS     ', 'ACS_S1       ',  &
    'ACS_S9       ', 'ACS_DIAMOND  ', 'ACS_CKBOARD  ', 'ACS_DEGREE   ',  &
    'ACS_PLMINUS  ', 'ACS_BULLET   ', 'ACS_LARROW   ', 'ACS_RARROW   ',  &
    'ACS_UARROW   ', 'ACS_DARROW   ', 'ACS_BOARD    ', 'ACS_LANTERN  ',  &
    'ACS_BLOCK    '  /)
  integer (C_LONG), dimension(ACSNUM) :: acs_values = (/  &
    acs_ulcorner , acs_urcorner , acs_llcorner, acs_lrcorner ,  &
    acs_ltee     , acs_rtee     , acs_ttee    , acs_btee     ,  &
    acs_hline    , acs_vline    , acs_plus    , acs_s1       ,  &
    acs_s9       , acs_diamond  , acs_ckboard , acs_degree   ,  &
    acs_plminus  , acs_bullet   , acs_larrow  , acs_rarrow   ,  &
    acs_uarrow   , acs_darrow   , acs_board   , acs_lantern  ,  &
    acs_block   /)
  integer :: i
  integer (C_INT) :: ierr, ix, iy, tmarg

  tmarg = (LINES - 22) / 2;

  ierr=attrset(curses_a_bold)
  ierr= mvaddstr(tmarg, (COLS - 23) / 2,  &
                 'Alternate Character Set' // C_NULL_CHAR)
  ierr=attrset(curses_a_normal)

  tmarg=tmarg+3

  do i=1, ACSNUM
    iy=mod(i-1,8) * 2 + tmarg
    ix=((i-1) / 8) * (COLS / 4) + (COLS / 8 - 7)
    ierr=move(iy, ix)
    ierr=addch(acs_values(i))
    ierr=mvaddstr(iy, ix+2, trim(acs_names(i)) // C_NULL_CHAR)
  end do
  ierr=mvaddstr(tmarg + 18, 3, 'Press any key to continue' // C_NULL_CHAR)
  ierr=getch()
end subroutine acstest

subroutine notavail(iwin)
  use curses
  use winsize
  type (C_PTR), intent(inout) :: iwin
  integer (C_INT) :: ierr

  ierr=beep()
  ierr=werase(iwin)
  ierr=box(iwin, acs_vline, acs_hline)
  ierr=cbreak()
  ierr=mvwaddstr(iwin, 2, 5,  &
        'Sorry, not yet implemented' // C_NULL_CHAR)
  call cont(iwin)
end subroutine notavail


subroutine cont(iwin)
  use curses
  type (C_PTR), intent(inout) :: iwin
  integer (C_INT) :: ierr
  ierr=mvwaddstr(iwin, 10, 1, ' Press any key to continue' // C_NULL_CHAR)
  ierr=wrefresh(iwin)
  ierr=raw()
  ierr=wgetch(iwin)
end subroutine cont

subroutine cont2()
  use curses
  integer (C_INT) :: ierr
  ierr=move(LINES-1,1)
  ierr=clrtoeol()
  ierr=mvaddstr(LINES-2, 1, ' Press any key to continue' // C_NULL_CHAR)
  ierr=refresh()
  ierr=raw()
  ierr=getch()
end subroutine cont2

subroutine display_menu(old_option, new_option)
  use curses
  use commands
  integer, intent(in) :: old_option, new_option
  integer (C_INT) :: lmarg, tmarg
  integer (C_INT) :: ierr
  integer :: i

  lmarg = (COLS - 14)/2
  tmarg = (LINES - (MAX_OPTIONS + 2))/2
  if (old_option == 0) then
    ierr=attrset(curses_a_bold)
    ierr=mvaddstr(tmarg-3, lmarg-5, 'PDCurses Test Program' // C_NULL_CHAR)
    ierr=attrset(curses_a_normal)
    do i=1, MAX_OPTIONS
      ierr=mvaddstr(tmarg+i, lmarg, trim(command(i)) // C_NULL_CHAR)
    end do
  else
    ierr=mvaddstr(tmarg + old_option, lmarg,  &
                  trim(command(old_option)) // C_NULL_CHAR)
  end if
  ierr=attrset(curses_a_reverse)
  ierr=mvaddstr(tmarg + new_option, lmarg,  &
                trim(command(new_option)) // C_NULL_CHAR)
  ierr=attrset(curses_a_normal)
  ierr=mvaddstr(tmarg + MAX_OPTIONS + 2, lmarg - 23,  &
         "Use Up and Down Arrows to select - Enter to run - Q to quit"  &
             // C_NULL_CHAR)
    ierr=refresh()
end subroutine display_menu

program testcurses
  use curses
  use commands
  type (C_PTR) :: iwin = C_NULL_PTR
  integer (C_INT) :: key
  integer :: istat, new_option=1, old_option=0
  call initTest(iwin, istat)
  if (istat /= 0) then
    write(*,'(a)') 'ERROR: initscr failed!'
    stop
  end if
! ierr=init_pair(1, curses_color_white, curses_color_blue)
! ierr=wbkgd(iwin,
  ierr=wbkgd(iwin, curses_a_reverse)
  ierr=erase()
  call display_menu(old_option, new_option)
  do
    ierr=noecho()
    ierr=raw()
    ierr=keypad(stdscr, curses_true)
    key=getch()
    if (key == 10 .or. key == 13 .or. key==curses_key_enter) then
      old_option = 0
      ierr=erase()
      ierr=refresh()
      if (new_option == 1) then
        call introtest(iwin, istat)
      else if (new_option == 2) then
        call resizetest(istat)
      else if (new_option == 3) then
        call scrolltest(iwin, istat)
      else if (new_option == 5) then
        call outputtest(iwin)
      else if (new_option == 6) then
        call acstest(iwin)
      else
        call notavail(iwin)
      end if
      ierr=erase()
      call display_menu(old_option, new_option)
    else if (key==curses_key_ppage .or. key==curses_key_home) then
      old_option = new_option
      new_option = 1
      call display_menu(old_option, new_option)
    else if (key==curses_key_npage .or. key==curses_key_end) then
      old_option = new_option
      new_option = MAX_OPTIONS
      call display_menu(old_option, new_option)
    else if (key==curses_key_up) then
      old_option = new_option
      if (new_option > 1) then
        new_option = new_option-1
      end if
      call display_menu(old_option, new_option)
    else if (key==curses_key_down) then
      old_option = new_option
      if (new_option < MAX_OPTIONS) then
        new_option = new_option+1
      end if
      call display_menu(old_option, new_option)
    else if (key==ichar('Q') .or. key==ichar('q')) then
      exit
    end if
  end do
  ierr=delwin(iwin)
  ierr=endwin()
end program testcurses
er (C_LONG), value, intent(in) :: attrs
    end function slk_attrset
!         int slk_attr_set(const attr_t attrs, short color_pair, void *opts);
    function slk_attr_set(attrs,color_pair,opts) result (ires) bind(C, name='slk_attr_set')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
      integer (C_SHORT), value :: color_pair
      type (C_PTR), value :: opts
    end function slk_attr_set
!         int slk_attroff(const chtype attrs);
    function slk_attroff(attrs) result (ires) bind(C, name='slk_attroff')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
    end function slk_attroff
!         int slk_attr_off(const attr_t attrs, void *opts);
    function slk_attr_off(attrs,opts) result (ires) bind(C, name='slk_attr_off')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_LONG), value, intent(in) :: attrs
      type (C_PTR), value :: opts
    end function slk_attr_off
!         int slk_color(short color_pair);
    function slk_color(color_pair) result (ires) bind(C, name='slk_color')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_SHORT), value :: color_pair
    end function slk_color
!         int PDC_mouse_in_slk(int y, int x);
    function PDC_mouse_in_slk(y,x) result (ires) bind(C, name='PDC_mouse_in_slk')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function PDC_mouse_in_slk
!         void PDC_slk_free(void);
    subroutine PDC_slk_free() bind(C, name='PDC_slk_free')
      use iso_c_binding
    end subroutine PDC_slk_free
!         void PDC_slk_initialize(void);
    subroutine PDC_slk_initialize() bind(C, name='PDC_slk_initialize')
      use iso_c_binding
    end subroutine PDC_slk_initialize
!         int baudrate(void);
    function baudrate() result (ires) bind(C, name='baudrate')
      use iso_c_binding
      integer (C_INT) :: ires
    end function baudrate
!         char erasechar(void);
    function erasechar() result (ich) bind(C, name='erasechar')
      use iso_c_binding
      character (C_CHAR) :: ich
    end function erasechar
!         bool has_ic(void);
    function has_ic() result (ires) bind(C, name='has_ic')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function has_ic
!         bool has_il(void);
    function has_il() result (ires) bind(C, name='has_il')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
    end function has_il
!         char killchar(void);
    function killchar() result (ich) bind(C, name='killchar')
      use iso_c_binding
      character (C_CHAR) :: ich
    end function killchar
!         char *longname(void);
!         chtype termattrs(void);
    function termattrs() result (ires) bind(C, name='termattrs')
      use iso_c_binding
      integer (C_LONG) :: ires
    end function termattrs
!         attr_t term_attrs(void);
    function term_attrs() result (ires) bind(C, name='term_attrs')
      use iso_c_binding
      integer (C_LONG) :: ires
    end function term_attrs
!         char *termname(void);
!         char wordchar(void);
    function wordchar() result (ich) bind(C, name='wordchar')
      use iso_c_binding
      character (C_CHAR) :: ich
    end function wordchar
!         int mvcur(int oldrow, int oldcol, int newrow, int newcol);
    function mvcur(oldrow,oldcol,newrow,newcol) result (ires) bind(C, name='mvcur')
      use iso_c_binding
      integer (C_INT) :: ires
      integer (C_INT), value :: oldrow
      integer (C_INT), value :: oldcol
      integer (C_INT), value :: newrow
      integer (C_INT), value :: newcol
    end function mvcur

!         int touchwin(WINDOW *win);
    function touchwin(win) result (ires) bind(C, name='touchwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function touchwin
!         int touchline(WINDOW *win, int start, int count);
    function touchline(win,start,cnt) result (ires) bind(C, name='touchline')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: start
      integer (C_INT), value :: cnt
    end function touchline
!         int untouchwin(WINDOW *win);
    function untouchwin(win) result (ires) bind(C, name='untouchwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function untouchwin
!         int wtouchln(WINDOW *win, int y, int n, int changed);
    function wtouchln(win,y,n,changed) result (ires) bind(C, name='wtouchln')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: n
      integer (C_INT), value :: changed
    end function wtouchln
!         bool is_linetouched(WINDOW *win, int line);
    function is_linetouched(win,line) result (ires) bind(C, name='is_linetouched')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: line
    end function is_linetouched
!         bool is_wintouched(WINDOW *win);
    function is_wintouched(win) result (ires) bind(C, name='is_wintouched')
      use iso_c_binding
      integer (C_SIGNED_CHAR) :: ires
      type (C_PTR), value :: win
    end function is_wintouched
!         WINDOW *newwin(int nlines, int ncols, int begy, int begx);
    function newwin(nlines,ncols,begy,begx) result (iwin) bind(C, name='newwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function newwin
!       WINDOW *derwin(WINDOW* orig, int nlines, int ncols, int begy, int begx);
    function derwin(orig,nlines,ncols,begy,begx) result (iwin) bind(C, name='derwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: orig
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function derwin
!       WINDOW *subwin(WINDOW* orig, int nlines, int ncols, int begy, int begx);
    function subwin(orig,nlines,ncols,begy,begx) result (iwin) bind(C, name='subwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: orig
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function subwin
!         WINDOW *dupwin(WINDOW *win);
    function dupwin(win) result (iwin) bind(C, name='dupwin')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: win
    end function dupwin
!         int delwin(WINDOW *win);
    function delwin(win) result (ires) bind(C, name='delwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
    end function delwin
!         int mvwin(WINDOW *win, int y, int x);
    function mvwin(win,y,x) result (ires) bind(C, name='mvwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: y
      integer (C_INT), value :: x
    end function mvwin
!         int mvderwin(WINDOW *win, int pary, int parx);
    function mvderwin(win,pary,parx) result (ires) bind(C, name='mvderwin')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: pary
      integer (C_INT), value :: parx
    end function mvderwin
!         int syncok(WINDOW *win, bool bf);
    function syncok(win, bf) result (ires) bind(C, name='syncok')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_SIGNED_CHAR), value :: bf
    end function syncok
!         void wsyncup(WINDOW *win);
    subroutine wsyncup(win) bind(C, name='wsyncup')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine wsyncup
!         void wcursyncup(WINDOW *win);
    subroutine wcursyncup(win) bind(C, name='wcursyncup')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine wcursyncup
!         void wsyncdown(WINDOW *win);
    subroutine wsyncdown(win) bind(C, name='wsyncdown')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine wsyncdown

!         WINDOW *resize_window(WINDOW *win, int nlines, int ncols);
    function resize_window(win,nlines,ncols) result (iwin) bind(C, name='resize_window')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: win
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
    end function resize_window
!         int wresize(WINDOW *win, int nlines, int ncols);
    function wresize(win,nlines,ncols) result (ires) bind(C, name='wresize')
      use iso_c_binding
      integer (C_INT) :: ires
      type (C_PTR), value :: win
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
    end function wresize
!         WINDOW *PDC_makelines(WINDOW *win);
    function PDC_makelines(win) result (iwin) bind(C, name='PDC_makelines')
      use iso_c_binding
      type (C_PTR) :: iwin
      type (C_PTR), value :: win
    end function PDC_makelines
!         WINDOW *PDC_makenew(int nlines, int ncols, int begy, int begx);
    function PDC_makenew(nlines,ncols,begy,begx) result (iwin) bind(C, name='PDC_makenew')
      use iso_c_binding
      type (C_PTR) :: iwin
      integer (C_INT), value :: nlines
      integer (C_INT), value :: ncols
      integer (C_INT), value :: begy
      integer (C_INT), value :: begx
    end function PDC_makenew
!         void PDC_sync(WINDOW *win);
    subroutine PDC_sync(win) bind(C, name='PDC_sync')
      use iso_c_binding
      type (C_PTR), value :: win
    end subroutine PDC_sync
  end interface
contains
!         int getch(void);
  function getch() result(ires)
    integer (C_INT) :: ires
    ires=wgetch(stdscr)
  end function getch
end module ncurses
