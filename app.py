import os, requests, streamlit as st, smtplib, time
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from typing import List, Dict
from Bio import Entrez
from langchain_openai import AzureChatOpenAI

# ░ ENV vars (set in Streamlit Cloud ▸ Secrets)
AZURE_ENDPOINT    = os.environ["AZURE_OPENAI_ENDPOINT"]
AZURE_KEY         = os.environ["AZURE_OPENAI_KEY"]
BREVO_SMTP_USER   = os.environ["BREVO_SMTP_USER"]
BREVO_SMTP_PASS   = os.environ["BREVO_SMTP_PASS"]
BREVO_SENDER_EMAIL= os.environ["BREVO_SENDER_EMAIL"]
BREVO_SENDER_NAME = os.environ["BREVO_SENDER_NAME"]

Entrez.email = "you@example.com"

llm = AzureChatOpenAI(
    azure_endpoint     = AZURE_ENDPOINT,
    api_key            = AZURE_KEY,
    azure_deployment   = "gpt-4o",
    openai_api_version = "2024-12-01-preview",
    temperature        = 0,
)

# ░ Helper functions
def fetch_pubmed(term:str,n:int,look:int|None)->List[Dict]:
    args=dict(db="pubmed",term=term,retmax=n)
    if look:
        args.update(mindate=(datetime.utcnow()-timedelta(days=look)).strftime("%Y/%m/%d"),
                    datetype="pdat")
    ids=Entrez.read(Entrez.esearch(**args))["IdList"]
    if not ids:return[]
    docs=Entrez.read(Entrez.esummary(db="pubmed",id=",".join(ids)))
    return[{"title":d["Title"][:250],"date":d.get("PubDate","N/A"),
            "url":f"https://pubmed.ncbi.nlm.nih.gov/{d['Id']}/"}for d in docs]

def fetch_trials(cond:str,kw:str,n:int)->List[Dict]:
    r=requests.get("https://clinicaltrials.gov/api/v2/studies",params={
        "query.cond":cond,"query.term":kw,"pageSize":n,"format":"json"},timeout=20)
    r.raise_for_status()
    out=[]
    for s in r.json().get("studies",[]):
        ps=s["protocolSection"]; ident=ps["identificationModule"]
        phase=ps.get("designModule",{}).get("phaseList",{}).get("phase",["N/A"])[0]
        out.append({"title":ident.get("officialTitle")or ident["briefTitle"],
                    "nct":ident["nctId"],"phase":phase,
                    "start":ps["statusModule"].get("startDateStruct",{}).get("startDate","N/A"),
                    "url":f"https://clinicaltrials.gov/study/{ident['nctId']}"})
    return out

def bullet(txt:str)->str:
    return llm.invoke(
        f"Summarize in ≤30 words for oncologists:\n\n{txt}\n\nBullet:"
    ).content.strip().lstrip("- ").strip()

def send_email(to:str,sub:str,body:str):
    msg=MIMEText(body); msg["Subject"]=sub
    msg["From"]=f"{BREVO_SENDER_NAME} <{BREVO_SENDER_EMAIL}>"; msg["To"]=to
    with smtplib.SMTP("smtp-relay.brevo.com",587) as s:
        s.starttls(); s.login(BREVO_SMTP_USER,BREVO_SMTP_PASS); s.send_message(msg)

# ░ Streamlit UI
st.title("🧬 Oncology Digest Generator")

topic    = st.text_input("Search term","ctDNA bladder cancer")
rec_name = st.text_input("Recipient name","Dr. Bukavina")
email    = st.text_input("Recipient email")
papers_n = st.slider("Max papers",1,100,8)
trials_n = st.slider("Max trials",1,50,5)
days     = {"7 days":7,"30 days":30,"All time":None}[st.selectbox(
            "PubMed window",("7 days","30 days","All time"))]

if st.button("Generate & Send"):

    progress = st.progress(0)
    with st.status("Generating digest…",expanded=True) as status:

        status.write("🔍 Fetching PubMed articles…")
        papers=fetch_pubmed(topic,papers_n,days); progress.progress(20)

        status.write("🔍 Fetching clinical trials…")
        trials=fetch_trials("bladder cancer",topic,trials_n); progress.progress(40)

        status.write("🧠 Summarizing with GPT…")
        for p in papers: p["summary"]=bullet(p["title"]); progress.progress(60)
        for t in trials: t["summary"]=bullet(t["title"]); progress.progress(80)

        today=datetime.now().strftime("%d %b %Y")
        body=["Dear "+rec_name+",","",
              f'Here is your digest on "{topic}":',"",
              "📚 Papers:"]
        body+= (sum([[f"• {p['summary']}",f"  ({p['date']}) {p['url']}"] for p in papers],[])
                or ["  – None found."])
        body+= ["","🧪 Trials:"]
        body+= (sum([[f"• {t['summary']} [{t['phase']}]",
                      f"  NCT {t['nct']} — Start {t['start']}",
                      f"  {t['url']}"] for t in trials],[])
                or ["  – None found."])
        body+= ["","Regards,","Oncology AI"]
        email_body="\n".join(body)
        subject=f"Oncology Digest – {topic} ({today})"
        progress.progress(90)

        if email:
            try:
                send_email(email,subject,email_body)
                status.write("✅ Email sent!")
            except Exception as e:
                status.error(f"Email failed: {e}")
        else:
            status.warning("No email entered; digest not sent.")

        progress.progress(100)
        status.write("🚀 Done!")

    st.text_area("Digest preview",email_body,height=400)
