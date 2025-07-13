import os, requests, streamlit as st, smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from typing import List, Dict
from Bio import Entrez
from langchain_openai import AzureChatOpenAI

# ── SECRETS (set in Streamlit Cloud) ───────────────────────────
AZURE_ENDPOINT       = os.environ["AZURE_OPENAI_ENDPOINT"]
AZURE_KEY            = os.environ["AZURE_OPENAI_KEY"]
BREVO_SMTP_USER      = os.environ["BREVO_SMTP_USER"]
BREVO_SMTP_PASS      = os.environ["BREVO_SMTP_PASS"]
BREVO_SENDER_EMAIL   = os.environ["BREVO_SENDER_EMAIL"]
BREVO_SENDER_NAME    = os.environ["BREVO_SENDER_NAME"]

Entrez.email = "you@example.com"  # PubMed contact

llm = AzureChatOpenAI(
    azure_endpoint     = AZURE_ENDPOINT,
    api_key            = AZURE_KEY,
    azure_deployment   = "gpt-4o",
    openai_api_version = "2024-12-01-preview",
    temperature        = 0,
)

# ── HELPERS ────────────────────────────────────────────────────
def fetch_pubmed(term: str, max_results: int, lookback: int | None) -> List[Dict]:
    args = dict(db="pubmed", term=term, retmax=max_results)
    if lookback:
        since = (datetime.utcnow() - timedelta(days=lookback)).strftime("%Y/%m/%d")
        args.update(mindate=since, datetype="pdat")
    ids = Entrez.read(Entrez.esearch(**args))["IdList"]
    if not ids: return []
    docs = Entrez.read(Entrez.esummary(db="pubmed", id=",".join(ids)))
    return [{"title": d["Title"][:250], "date": d.get("PubDate","N/A"),
             "url":f"https://pubmed.ncbi.nlm.nih.gov/{d['Id']}/"} for d in docs]

def fetch_trials(cond: str, kw: str, n: int) -> List[Dict]:
    r = requests.get("https://clinicaltrials.gov/api/v2/studies", params={
        "query.cond": cond, "query.term": kw,
        "pageSize": n, "format":"json"}, timeout=20)
    r.raise_for_status()
    items = []
    for s in r.json().get("studies", []):
        ps = s["protocolSection"]
        items.append({
            "title":  ps["identificationModule"].get("officialTitle") or
                      ps["identificationModule"]["briefTitle"],
            "nct":    ps["identificationModule"]["nctId"],
            "phase":  ps.get("designModule",{}).get("phaseList",{}).get("phase",["N/A"])[0],
            "start":  ps["statusModule"].get("startDateStruct",{}).get("startDate","N/A"),
            "url":    f"https://clinicaltrials.gov/study/{ps['identificationModule']['nctId']}"
        })
    return items

def bullet(txt:str)->str:
    summary=llm.invoke(
        f"Summarize in ≤30 words for oncologists:\n\n{txt}\n\nBullet:"
    ).content.strip()
    return summary.lstrip("- ").strip()  # remove leading dash if GPT adds one

def send_email(to_mail:str, subject:str, body:str):
    msg=MIMEText(body)
    msg["Subject"]=subject
    msg["From"]=f"{BREVO_SENDER_NAME} <{BREVO_SENDER_EMAIL}>"
    msg["To"]=to_mail
    with smtplib.SMTP("smtp-relay.brevo.com",587) as s:
        s.starttls(); s.login(BREVO_SMTP_USER,BREVO_SMTP_PASS); s.send_message(msg)

# ── STREAMLIT UI ───────────────────────────────────────────────
st.title("🧬 Oncology Digest Generator")

topic      = st.text_input("Search term", "ctDNA bladder cancer")
rec_name   = st.text_input("Recipient name", "Dr. Bukavina")
recipient  = st.text_input("Recipient email")
max_papers = st.slider("Max papers", 1, 100, 8)
max_trials = st.slider("Max trials", 1, 50, 5)
lookback   = st.selectbox("PubMed window", ("7 days","30 days","All time"))
days       = {"7 days":7,"30 days":30,"All time":None}[lookback]

if st.button("Generate & Send"):
    with st.spinner("Collecting literature…"):
        papers = fetch_pubmed(topic, max_papers, days)
        trials = fetch_trials("bladder cancer", topic, max_trials)

    for p in papers:  p["summary"]=bullet(p["title"])
    for t in trials:  t["summary"]=bullet(t["title"])

    today = datetime.now().strftime("%d %b %Y")
    body_lines=[
        f"Dear {rec_name},",
        "",
        f'Here is your digest on "{topic}":',
        "",
        "📚 Papers:"
    ]+(
        sum([[f"• {p['summary']}",f"  ({p['date']}) {p['url']}"] for p in papers],[])
        or ["  – None found."]
    )+[
        "",
        "🧪 Trials:"
    ]+(
        sum([[f"• {t['summary']} [{t['phase']}]",
              f"  NCT {t['nct']} — Start {t['start']}",
              f"  {t['url']}"] for t in trials],[])
        or ["  – None found."]
    )+["","Regards,","Oncology AI"]

    email_body="\n".join(body_lines)
    email_subject=f"Oncology Digest – {topic} ({today})"

    st.text_area("Digest preview", email_body, height=400)

    if recipient:
        try:
            send_email(recipient, email_subject, email_body)
            st.success(f"Email sent to {recipient}")
        except Exception as e:
            st.error(f"Email failed: {e}")
    else:
        st.warning("Enter a recipient email to send.")
