"""
Here we provide a simple Python script that can read the arXiv and email to you
the filtered abstracts according to some keywords that you determine.
"""
#!/usr/local/bin/python

# Import stuff
import requests
from bs4 import BeautifulSoup
import smtplib
from email.mime.text import MIMEText

# Define keywords. You can add as many as you want
keywords = ['keyword1', 'keyword2', 'keyword3']

# Define website to be read and read it correctly
link = "https://arxiv.org/list/cond-mat/new"
page = requests.get(link)
soup = BeautifulSoup( page.content, 'html.parser')

# Extract the information needed and put it in format to be checked
titles = soup.find_all('div', {'class' : 'list-title mathjax'})
abstracts = soup.find_all('p', {'class' : 'mathjax'})
authors = soup.find_all('div', {'class' : 'list-authors'})
refs = soup.find_all('a', {'title' : 'Abstract'})
lines_titles = [title.get_text() for title in titles]
lines_abstracts = [abstract.get_text() for abstract in abstracts]
lines_authors = [author.get_text() for author in authors]
lines_refs = [ref.get_text() for ref in refs]

# Write filtered papers on a file in your local computer.
# In this case /PathToWhateverYouWant/arxiv_summary.txt
filetosend = open('/PathToWhateverYouWant/arxiv_summary.txt','w')

for i in range(len(lines_abstracts)):
    if any(word in lines_abstracts[i] for word in keywords):
        filetosend.write(lines_titles[i].encode('ascii', 'ignore').decode('ascii'))
        filetosend.write(lines_authors[i].encode('ascii', 'ignore').decode('ascii'))
        filetosend.write('Abstract: ')
        filetosend.write(lines_abstracts[i].encode('ascii', 'ignore').decode('ascii'))
        filetosend.write('https://arxiv.org/')
        filetosend.write(lines_refs[i].encode('ascii', 'ignore'). \
                         decode('ascii').replace('arXiv:', 'abs/'))
        filetosend.write('\n**********\n')

filetosend.close()

# Prepare to send email. Read file and generate plain text.
fp = open(filetosend, 'rb')
msg = MIMEText(fp.read())
fp.close()

# Define Subject, From, To and send message
# Use the same email address to send and receive
msg['Subject'] = 'arXiv'
msg['From'] = 'email@youruni'
msg['To'] = 'email@youruni'

# Include your SMTP server address
s = smtplib.SMTP('smtpserver')
s.sendmail('email@youruni', 'email@youruni', msg.as_string())
s.quit()
"""
This script can be scheduled to be executed automatically using crontab. For instance, if you want to get your arXiv email every weekday at 12:00 use the command


    crontab -e

and include this line in the file that will appear


    0 12 * * 1-5 /PathToTheAboveScript.py

"""
