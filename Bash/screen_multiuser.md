screen which has multiuser support.

First, create a new session:

*screen -d -m -S multisession*

Attach to it:

*screen -r multisession*

Turn on multiuser support:

Press *Ctrl-a* and type (NOTE: Ctrl+a is needed just before each single command, i.e. twice here)

*:multiuser on*
*:acladd USER* ← use username of user you want to give access to your screen

Now, *Ctrl-a d* and list the sessions:

*$ screen -ls*
There is a screen on:
    *4791.multisession   (Multi, detached)*

You now have a multiuser screen session. Give the name multisession to acl'd user, so he can attach to it:

*screen -x youruser/multisession*

And that's it.

The only drawback is that screen must run as suid root. But as far as I know is the default, normal situation.

Another option is to do *screen -S $screen_id -X multiuser on, screen -S $screen_id -X acladd authorized_user*
