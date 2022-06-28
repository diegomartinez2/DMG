
subroutine c_f_string_chars(c_string, f_string)
    ! Helper function
    use iso_c_binding
    implicit none
    character(len=1,kind=c_char), intent(in) :: c_string(*)
    character(len=*), intent(out) :: f_string
    integer :: i
    i=1
    do while(c_string(i)/=c_null_char .and. i<=len(f_string))
        f_string(i:i) = c_string(i)
        i=i+1
    end do
    if (i<=len(f_string)) f_string(i:) = ' '
end subroutine c_f_string_chars
!
!#######################################################
!
module say_hello_m
    use iso_c_binding
    use gtk
    use gtk_sup
    use g
    implicit none

    type say_hello_t
        private
        type(c_ptr) :: window
        type(c_ptr) :: te_name
    end type

    private
    public :: say_hello_t, say_hello_new, say_hello_show

contains

    subroutine say_hello_new(this)
        type(say_hello_t), target :: this

        type(c_ptr)    :: builder
        type(c_ptr)    :: error
        integer(c_int) :: guint


        ! load GUI into builder
        builder = gtk_builder_new()
        error = c_null_ptr
        guint = gtk_builder_add_from_file(builder, "say_hello_gtk2.glade"//c_null_char, error)
        if (guint == 0) then
            print *, "Could not open say_hello.glade"
            stop "Program terminated"
        end if

        ! get references to GUI elements
        this%window  = gtk_builder_get_object(builder,"window"//c_null_char)
        this%te_name = gtk_builder_get_object(builder,"te_name"//c_null_char)

        ! connect signal handlers
        call gtk_builder_connect_signals_full(builder, c_funloc(connect_signals), c_loc(this))

        ! free memory
        call g_object_unref(builder)

    end subroutine

    !-----------------------------------------
    subroutine connect_signals (builder, object, signal_name, handler_name, connect_object, flags, c_this) bind(c)
        use iso_c_binding, only: c_ptr, c_char, c_int
        type(c_ptr), value                     :: builder        !a GtkBuilder
        type(c_ptr), value                     :: object         !object to connect a signal to
        character(kind=c_char), dimension(*)   :: signal_name    !name of the signal
        character(kind=c_char), dimension(*), target   :: handler_name   !name of the handler
        type(c_ptr), value                     :: connect_object !a GObject, if non-NULL, use g_signal_connect_object()
        integer(c_int), value                  :: flags          !GConnectFlags to use
        type(c_ptr), value                     :: c_this         !user data

        character(len=32)                      :: h_name


        call c_f_string_chars(handler_name, h_name)
        print *, "Connect signal for = ", h_name

        select case (h_name)
        case ("sh_bt_hello_clicked")
            call g_signal_connect (object, signal_name, c_funloc(sh_bt_hello_clicked), c_this)
        case ("sh_bt_quit_clicked")
            call g_signal_connect (object, signal_name, c_funloc(sh_bt_quit_clicked), c_this)
        case ("sh_quit_window")
            call g_signal_connect (object, signal_name, c_funloc(sh_quit_window), c_this)
        case default
            print *, "Unknow handler = "//h_name
            stop "Program terminated"
        end select

    end subroutine

    !-----------------------------------------
    subroutine say_hello_show(this)
        type(say_hello_t), intent(in) :: this

        call gtk_widget_show(this%window)
    end subroutine

    !-----------------------------------------
    function sh_bt_hello_clicked(widget, gdata) result(ret) bind(c)
        use gtk_hl_dialog

        type(c_ptr), value :: widget, gdata
        integer(c_int) :: ret
        integer(c_int) :: resp

        type(say_hello_t), pointer :: this
        integer(c_int16_t) :: elen
        type(c_ptr) :: cptxt
        character(kind=c_char,len=100), pointer :: cname
        character(kind=c_char,len=30)  :: name
        character(len=80) :: msg

        call c_f_pointer(gdata, this)

        print *, "Say hello was clicked"

        elen = gtk_entry_get_text_length(this%te_name)
        print *, 'elen=', elen
        if (elen == 0) then
            msg = "Write your name first, please."
        else
            cptxt = gtk_entry_get_text(this%te_name)
            call c_f_string(cptxt,name)
            msg  = "Hello "//trim(name)//"!"
        end if

        resp = hl_gtk_message_dialog_show([msg], GTK_BUTTONS_OK, "Message"//c_null_char, parent=this%window)

        ret = FALSE
    end function

    !-----------------------------------------
    function sh_bt_quit_clicked(widget, gdata) result(ret) bind(c)
        type(c_ptr) :: widget, gdata
        integer(c_int) :: ret

        print *, "Quit was clicked"
        call gtk_main_quit()
        ret = FALSE
    end function

    !-----------------------------------------
    function sh_quit_window(widget, gdata) result(ret) bind(c)
        type(c_ptr) :: widget, gdata
        integer(c_int) :: ret

        print *, "Quit window"
        call gtk_main_quit()
        ret = FALSE
    end function

end module

program say_hello_main
    use gtk, only: gtk_init, gtk_main
    use say_hello_m
    implicit none

    type(say_hello_t) :: say_hello_win

    call gtk_init()

    call say_hello_new(say_hello_win)
    call say_hello_show(say_hello_win)

    call gtk_main()
end program say_hello_main
